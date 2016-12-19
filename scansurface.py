"""

Copyright (C) 2016 Jan Jaeken <jan.jaeken@ugent.be>

This file is part of Christoffel.

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

#Read special directions if present
try:
    directions = map(float, config.get('SCAN', 'directions').split())
    directions = np.reshape(directions, (len(directions)/3, 3))
except:
    directions = []

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

#Prepare the output files
SS_outfile = open('slow_secondary.dat', 'w')
SS_outfile.write('# Slow secondary\n')
SS_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
SS_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
SS_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
SS_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

FS_outfile = open('fast_secondary.dat', 'w')
FS_outfile.write('# Fast secondary\n')
FS_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
FS_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
FS_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
FS_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

P_outfile = open('primary.dat', 'w')
P_outfile.write('# Primary\n')
P_outfile.write('# Theta (rad) # Phi (rad) # Cube x y z ')
P_outfile.write('# Phase (km/s) # Phase (rel) # Phase pol x y z ')
P_outfile.write('# Group (km/s) # Group (rel) # Group x y z (km/s) ')
P_outfile.write('# Power flow angle (deg) # Enhancement factor\n')

all_out = [SS_outfile, FS_outfile, P_outfile]

#Prepare min/max stuff
maxvel = [1e-8, 1e-8, 1e-8]
minvel = [1e8, 1e8, 1e8]
maxdir = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
mindir = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
maxangle = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
minangle = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]

#Get the isotropic velocities
iso_vel = chris.get_isotropic()

print 'Calculating phase velocities in user-defined directions.'
#Loop over user-defined directions
direction_file = open('directions.dat', 'w')
for q in directions:
    chris.set_direction_cartesian(q)
    velocity = chris.get_phase_velocity()

    direction_file.write('Direction: ' + str(q[0]) + ' ' + str(q[1]) + ' ' + str(q[2]) + '\n')
    direction_file.write(' Primary: ' + str(velocity[2]) + ' km/s\n')
    direction_file.write(' Secondary 1: ' + str(velocity[1]) + ' km/s\n')
    direction_file.write(' Secondary 2: ' + str(velocity[0]) + ' km/s\n')
    direction_file.write(' Secondary avg: ' + str(0.5*(velocity[0]+velocity[1])) + ' km/s\n\n')
direction_file.close()
print 'The directions.dat file has been successfully written.'

print 'Looping over surface.'
#Loop over surface
for theta in xrange(num_theta):
    print str(theta+1) + '/' + str(num_theta)
    num_phi = total_num_phi
    #This commented bit scales num_phi by sin(theta) to get a constant sampling density
    #Sadly, GNUPlot doesn't play nice with grids that aren't rectangular
    #Still, if you want to cut your runtime down significantly, use the next two lines
    #sin_theta = np.sin(0.5*np.pi*theta / (num_theta - 1))
    #num_phi = int(np.ceil(sin_theta*total_num_phi)) + 1
    for phi in xrange(num_phi):
        chris.set_direction_spherical(0.5*np.pi*theta/(num_theta-1), 2.0*np.pi*phi/ (num_phi-1))
        q = chris.direction
        cubepos = q/np.linalg.norm(q, ord=np.inf)

        #Calculate everything and store in variables
        enhancefac = chris.get_enhancement()
        velocity = chris.get_phase_velocity()
        polarization = chris.get_eigenvec()
        group_vel = chris.get_group_velocity()
        group_abs = chris.get_group_abs()
        pf_angle = chris.get_powerflow()

        #Check for fastest and slowest velocities
        for i in xrange(3):
            if velocity[i] > maxvel[i]:
                maxvel[i] = velocity[i]
                maxdir[i] = q
                maxangle[i] = [chris.theta, chris.phi]
            if velocity[i] < minvel[i]:
                minvel[i] = velocity[i]
                mindir[i] = q
                minangle[i] = [chris.theta, chris.phi]

            #And start writing data
            all_out[i].write('{0:.6f} {1:.6f}  {2: 8.6f} {3: 8.6f} {4: 8.6f}  '.format(chris.theta, chris.phi, *cubepos))
            all_out[i].write('{0: 10.6f}  {1: 8.6f}  {2: 8.6f} {3: 8.6f} {4: 8.6f}  '.format(velocity[i], 100*(velocity[i]/iso_vel[i]-1.0), *(polarization[i])))
            all_out[i].write('{0: 10.6f}  {1: 8.6f}  {2: 10.6f} {3: 10.6f} {4: 10.6f}  '.format(group_abs[i], 100*(group_abs[i]/iso_vel[i]-1.0), *(group_vel[i])))
            all_out[i].write('{0: 8.4f}  {1:.7g}\n'.format(180*pf_angle[i]/np.pi, enhancefac[i]))

    SS_outfile.write('\n')
    FS_outfile.write('\n')
    P_outfile.write('\n')
SS_outfile.close()
FS_outfile.close()
P_outfile.close()

print 'Calculating anisotropy.'

#Rounding gives nicer output
maxdir = np.around(maxdir, 10)
mindir = np.around(mindir, 10)

#Write anisotropy data
anisotropy_file = open('anisotropy.dat', 'w')
anisotropy_file.write('Max Primary: ' + str(maxvel[2]) + ' km/s')
anisotropy_file.write(' at direction (' + str(maxangle[2][0]) + ', ' + str(maxangle[2][1]) + ')')
anisotropy_file.write(' - (' + str(maxdir[2][0]) + ', ' + str(maxdir[2][1]) + ', ' + str(maxdir[2][2]) + ')\n')
anisotropy_file.write('Min Primary: ' + str(minvel[2]) + ' km/s')
anisotropy_file.write(' at direction (' + str(minangle[2][0]) + ', ' + str(minangle[2][1]) + ')')
anisotropy_file.write(' - (' + str(mindir[2][0]) + ', ' + str(mindir[2][1]) + ', ' + str(mindir[2][2]) + ')\n\n')

anisotropy_file.write('Max Secondary 1: ' + str(maxvel[1]) + ' km/s')
anisotropy_file.write(' at direction (' + str(maxangle[1][0]) + ', ' + str(maxangle[1][1]) + ')')
anisotropy_file.write(' - (' + str(maxdir[1][0]) + ', ' + str(maxdir[1][1]) + ', ' + str(maxdir[1][2]) + ')\n')
anisotropy_file.write('Min Secondary 1: ' + str(minvel[1]) + ' km/s')
anisotropy_file.write(' at direction (' + str(minangle[1][0]) + ', ' + str(minangle[1][1]) + ')')
anisotropy_file.write(' - (' + str(mindir[1][0]) + ', ' + str(mindir[1][1]) + ', ' + str(mindir[1][2]) + ')\n\n')

anisotropy_file.write('Max Secondary 2: ' + str(maxvel[0]) + ' km/s')
anisotropy_file.write(' at direction (' + str(maxangle[0][0]) + ', ' + str(maxangle[0][1]) + ')')
anisotropy_file.write(' - (' + str(maxdir[0][0]) + ', ' + str(maxdir[0][1]) + ', ' + str(maxdir[0][2]) + ')\n')
anisotropy_file.write('Min Secondary 2: ' + str(minvel[0]) + ' km/s')
anisotropy_file.write(' at direction (' + str(minangle[0][0]) + ', ' + str(minangle[0][1]) + ')')
anisotropy_file.write(' - (' + str(mindir[0][0]) + ', ' + str(mindir[0][1]) + ', ' + str(mindir[0][2]) + ')\n\n')

anisotropy_file.write('Primary isotropic velocity: ' + str(chris.iso_P) + ' km/s\n')
anisotropy_file.write('Secondary isotropic velocity: ' + str(chris.iso_S) + ' km/s\n')
anisotropy_file.write('Primary Anisotropy: ' + str(100*(1.0 - minvel[2]/maxvel[2])) + '\n')
anisotropy_file.write('Secondary Anisotropy: ' + str(100*(1.0 - minvel[0]/maxvel[1])) + '\n')

anisotropy_file.close()

print 'All done!'
