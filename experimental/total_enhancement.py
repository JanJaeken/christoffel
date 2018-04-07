"""

Copyright (C) 2018 Jan Jaeken <jan.jaeken@gmail.com>

This file is part of Christoffel.

This file is an attempt at calculating the acoustic energy carried away from
a point source in different directions. It is untested and is not guaranteed
to give reliable results. Please do your own testing if you wish to use it.

"Solving the Christoffel equation: Phase and group velocities"
Computer Physics Communications, 10.1016/j.cpc.2016.06.014

Christoffel is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Christoffel is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Christoffel.  If not, see <http://www.gnu.org/licenses/>.

"""

print 'WARNING! This feature has not been tested. Do not trust these results.'

import christoffel
import numpy as np
import ConfigParser

print 'Reading data and settings from the sound.in file.'
config = ConfigParser.ConfigParser()
config.read('sound.in')

# An error in reading the tensor or density should crash the script.
# Can't do anything without a stiffness tensor.
stiffness_tensor = map(float, config.get('SCAN', 'stiffness').split())
stiffness_tensor = np.reshape(stiffness_tensor, (6, 6))

density = config.getfloat('SCAN', 'density') #kg/m^3

# Creation of the central Christoffel object
chris = christoffel.Christoffel(stiffness_tensor, density)

#Read in rotation information if present
try:
    zdir = map(float, config.get('SCAN', 'zdir').split())
except:
    zdir = None
try:
    xdir = map(float, config.get('SCAN', 'xdir').split())
except:
    xdir = None

#Read grid density for the scan
try:
    num_theta = config.getint('SCAN', 'numtheta')
    total_num_phi = config.getint('SCAN', 'numphi')
except:
    num_theta = 0
    total_num_phi = 0

print 'Data read and Christoffel object created.\n'

#Dump the data that has been read
statusfile = open('sound.out', 'w')

statusfile.write('Density: {0:.2f} kg/m^3\n\n'.format(chris.density))
statusfile.write('Stiffness tensor in GPa:\n')
for i in xrange(6):
    statusfile.write('{0:8.2f} {1:8.2f} {2:8.2f} {3:8.2f} {4:8.2f} {5:8.2f}\n'.format(*(chris.stiffness2D[i])))

chris.rotate_tensor(x_dir=xdir, z_dir=zdir)

statusfile.write('\nRotated stiffness tensor in GPa:\n')
for i in xrange(6):
    statusfile.write('{0:8.2f} {1:8.2f} {2:8.2f} {3:8.2f} {4:8.2f} {5:8.2f}\n'.format(*(chris.stiffness2D[i])))

statusfile.write('\nBulk modulus: {0:.2f} GPa\n'.format(chris.bulk))
statusfile.write('Shear modulus: {0:.2f} GPa\n'.format(chris.shear))
statusfile.write('Isotropic primary velocity: {0:.2f} km/s\n'.format(chris.iso_P))
statusfile.write('Isotropic secondary velocity: {0:.2f} km/s\n'.format(chris.iso_S))

statusfile.close()

#Prepare full grid
total_enhancement = np.zeros((3, num_theta + 1, total_num_phi + 1))
group_theta_index = np.empty(3)
group_phi_index = np.empty(3)

print 'Looping over surface.'
#Loop over surface
for theta in xrange(num_theta):
    print str(theta+1) + '/' + str(num_theta)
    sin_theta = np.sin(0.5*np.pi*theta/(num_theta-1))
    num_phi = max(int(np.ceil(sin_theta*total_num_phi)), 2)
    for phi in xrange(num_phi):
        chris.set_direction_spherical(0.5*np.pi*theta/(num_theta-1), 2.0*np.pi*phi/(num_phi-1))

        #Calculate everything and store in variables
        enhancefac = chris.get_enhancement()
        group_theta = chris.get_group_theta()
        group_phi = chris.get_group_phi()
        
        #Sum it all in proper bins
        for i in xrange(3):
            group_theta_index = int((num_theta-1) * 2.0*group_theta[i] / np.pi)
            if group_theta_index > num_theta-1:
                group_theta_index = num_theta-1 - group_theta_index % num_theta
            
            sin_theta_2 = np.sin(0.5*np.pi*group_theta_index/(num_theta-1))
            num_phi_2 = max(int(np.ceil(sin_theta_2*total_num_phi)), 2)
            group_phi_index = int((num_phi_2-1) * 0.5*group_phi[i] / np.pi) % num_phi_2
            total_enhancement[i][group_theta_index][group_phi_index] += enhancefac[i]


print 'Writing data.'
            
total_out = open('total_enhancement.dat', 'w')
total_out.write('# Theta (rad) # Phi (rad) # Cube x y z ')
total_out.write('# Enhancement factor S1 # Enhancement factor S2 # Enhancement factor P\n')

for theta in xrange(num_theta):
    print str(theta+1) + '/' + str(num_theta)
    theta_rad = 0.5*np.pi*theta/(num_theta-1)
    sin_theta = np.sin(theta_rad)
    num_phi = max(int(np.ceil(sin_theta*total_num_phi)), 2)
    for phi in xrange(num_phi):
        phi_rad = 2.0*np.pi*phi/(num_phi-1)
        chris.set_direction_spherical(theta_rad, phi_rad)
        q = chris.direction
        cubepos = q/np.linalg.norm(q, ord=np.inf)
        
        total_out.write('{0:.6f} {1:.6f}  {2: 8.6f} {3: 8.6f} {4: 8.6f}  '.format(theta_rad, phi_rad, *cubepos))
        for i in xrange(3):
            total_out.write('{0: 8.6f} '.format(total_enhancement[i][theta][phi]))
        total_out.write('\n')
    total_out.write('\n')

total_out.close()

print 'All done!'
