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

import numpy as np

# Definition of Voigt notation
VOIGT = {00: 0, 11: 1, 22: 2, 12: 3, 21: 3, 02: 4, 20: 4, 01: 5, 10: 5}
idmat = np.identity(3)

class Christoffel:
    """
    Contains all information about the material, such as
    density and stiffness tensor. Given a reciprocal vector
    (sound wave direction), it can produce phase and group
    velocities and associated enhancement factors.

    After initialization, set a wave vector direction with
    set_direction or set_direction_spherical, after which any and all
    information can be gained from the get_* functions. All calculations
    will be done on the fly on a need-to-know basis.

    Keyword arguments:
    stiffness -- 6x6 stiffness tensor in GPa
    density -- density of the material in kg/m^3
    """

    def __init__(self, stiffness, density):
        self.bulk = get_bulk(stiffness)
        self.shear = get_shear(stiffness)
        self.iso_P, self.iso_S = isotropic_velocities(self.bulk, self.shear, density)

        stiffness = 0.5 * ( stiffness + stiffness.T)
        self.stiffness = np.array(de_voigt(stiffness))
        self.stiffness *= 1000.0/density
        self.density = density

        self.hessian_mat = hessian_christoffelmat(self.stiffness)

        self.clear_direction()

    def clear_direction(self):
        """Clear all direction-dependent data"""
        self.direction = None

        self.theta = None
        self.phi = None

        self.christoffel = None
        self._grad_mat = None

        self._eigenval = None
        self._eigenvec = None
        self._grad_eig_val = None
        self._hessian_eig = None

        self._phase_velocity = None
        self._group_velocity = None
        self._group_abs = None
        self._group_dir = None
        self.group_theta = None
        self.group_phi = None
        self._powflow_angle = None
        self._cos_pf_angle = None
        self._enhancement = None

    def rotate_tensor(self, rot_mat=None, x_dir=None, z_dir=None):
        """
        Apply rotation defined by rot_mat to the rank-4 tensor.
        If no rot_mat is given, rotate the tensor to align
        the z-axis with z_dir and the x-axis with x_dir if provided.
        """
        self.clear_direction()
        if rot_mat is None:
            rot = idmat
            if z_dir is not None and x_dir is None:
                z_dir = z_dir / norm(z_dir)
                rot = get_rot_mat(z_dir, [0.0, 0.0, 1.0])
            if x_dir is not None and z_dir is None:
                x_dir = x_dir / norm(x_dir)
                rot = get_rot_mat(x_dir, [1.0, 0.0, 0.0])
            if x_dir is not None and z_dir is not None:
                x_dir = x_dir / norm(x_dir)
                z_dir = z_dir / norm(z_dir)

                x_dir -= np.dot(x_dir, z_dir) * z_dir # Gram-Schmidt
                x_dir = x_dir / norm(x_dir)
                y_dir = np.cross(z_dir, x_dir)

                rot = np.array([x_dir, y_dir, z_dir])
        else:
            rot = rot_mat
        for i in xrange(4):
            self.stiffness = np.tensordot(rot, self.stiffness, (1, i))

        self.hessian_mat = hessian_christoffelmat(self.stiffness)

    def set_direction_cartesian(self, direction):
        """
        Define a wave vector in cartesian coordinates.
        It is always explicitly normalized to lie on the unit sphere.
        """
        self.clear_direction()

        self.direction = direction / norm(direction)
        q = self.direction

        x = q[0]
        y = q[1]
        z = q[2]
        if z >= 1.0 or z <= -1.0:
            if z > 0.0:
                self.theta = 0.0
            else:
                self.theta = np.pi
            self.phi = 0.0
        else:
            self.theta = np.arccos(z)
            sin_theta = np.sqrt(1 - z**2)

            cos_phi = x/sin_theta

            self.phi = np.arccos(cos_phi)
            if y < 0.0:
                self.phi = 2.0*np.pi - self.phi

        self.christoffel = np.dot(q, np.dot(q, self.stiffness))

    def set_direction_spherical(self, theta, phi):
        """
        Define a wave vector in spherical coordinates (rad).
        Theta is the polar angle, phi the azimuth.
        x = cos(phi) * sin(theta)
        y = sin(phi) * sin(theta)
        z = cos(theta)
        """
        self.clear_direction()

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        x = cos_phi * sin_theta
        y = sin_phi * sin_theta
        z = cos_theta
        q = np.array([x, y, z])

        self.theta = theta
        self.phi = phi
        self.direction = q

        self.christoffel = np.dot(q, np.dot(q, self.stiffness))

    def set_direction_random(self):
        """
        Generates a random wave vector direction.
        The distribution is uniform across the unit sphere.
        """
        self.clear_direction()

        cos_theta = np.random.ranf()
        phi = 2.0 * np.pi * np.random.ranf()
        sin_theta = np.sqrt(1 - cos_theta**2)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        q = np.array([cos_phi*sin_theta, sin_phi*sin_theta, cos_theta])
        self.direction = q
        self.phi = phi
        self.theta = np.arccos(cos_theta)

        self.christoffel = np.dot(q, np.dot(q, self.stiffness))


    def get_bulk(self):
        return self.bulk

    def get_shear(self):
        return self.shear

    def get_isotropic(self):
        """
        Returns sound velocities as if the material was isotropic.
        """
        return np.array([self.iso_S, self.iso_S, self.iso_P])

    def get_isotropic_P(self):
        return self.iso_P

    def get_isotropic_S(self):
        return self.iso_S

    def get_stiffness(self):
        return self.stiffness

    def get_density(self):
        return self.density

    def get_direction(self):
        return self.direction

    def get_direction_spherical(self):
        return np.array([self.theta, self.phi])

    def get_christoffel_matrix(self):
        return self.christoffel

    def get_grad_mat(self):
        """
        Returns the gradient of the Christoffel matrix.
        d/dx_n M_ij = sum_k q_k * ( C_inkj + C_iknj )
        gradmat[n][i][j] =  d/dx_n M_ij (note the indices)
        """
        if self._grad_mat is None:
            self.set_grad_mat()
        return self._grad_mat

    def get_eigenval(self):
        """
        Returns the eigenvalues of the Christoffel matrix, sorted low to high.
        """
        if self._eigenval is None:
            self.set_phase_velocity()
        return self._eigenval

    def get_eigenvec(self):
        """
        Returns the eigenvectors of the Christoffel matrix,
        sorted from low to high eigenvalue.
        """
        if self._eigenvec is None:
            self.set_phase_velocity()
        return self._eigenvec

    def get_grad_eigenval(self):
        """Returns the gradient of the eigenvalues."""
        if self._grad_eig_val is None:
            self.set_group_velocity()
        return self._grad_eig_val

    def get_phase_velocity(self):
        if self._phase_velocity is None:
            self.set_phase_velocity()
        return self._phase_velocity

    def get_relative_phase_velocity(self):
        """Returns phase velocity / isotropic velocity."""
        if self._phase_velocity is None:
            self.set_phase_velocity()
        return self._phase_velocity / self.get_isotropic()

    def get_group_velocity(self):
        if self._group_velocity is None:
            self.set_group_velocity()
        return self._group_velocity

    def get_group_abs(self):
        if self._group_abs is None:
            self.set_group_velocity()
        return self._group_abs

    def get_relative_group_velocity(self):
        """Returns group velocity / isotropic velocity."""
        if self._group_abs is None:
            self.set_group_velocity()
        return self._group_abs / self.get_isotropic()

    def get_group_dir(self):
        if self._group_dir is None:
            self.set_group_velocity()
        return self._group_dir

    def get_group_theta(self):
        if self.group_theta is None:
            self.set_group_velocity()
        return self.group_theta

    def get_group_phi(self):
        if self.group_phi is None:
            self.set_group_velocity()
        return self.group_phi

    def get_powerflow(self):
        if self._powflow_angle is None:
            self.set_group_velocity()
        return self._powflow_angle

    def get_cos_powerflow(self):
        if self._cos_pf_angle is None:
            self.set_group_velocity()
        return self._cos_pf_angle

    def get_hessian_mat(self):
        """
        Returns the hessian of the Christoffel matrix.
        hessmat[i][j][k][l] = d^2 M_kl / dx_i dx_j  (note the indices).
        """
        return self.hessian_mat
        
    def get_hessian_eig(self):
        """
        Returns the hessian of the eigenvalues of the Christoffel matrix.
        Hessian[n][i][j] = d^2 lambda_n / dx_i dx_j
        """
        if self._hessian_eig is None:
            self.set_hessian_eig()
        return self._hessian_eig

    def get_enhancement(self, approx=False, num_steps=8, delta=1e-5):
        if self._enhancement is None:
            if approx is False:
                self.set_enhancement()
            else:
                self.set_enhancement_approx(num_steps, delta)
        return self._enhancement


    def set_phase_velocity(self):
        """
        Determine eigenvalues, eigenvectors of the Christoffel matrix,
        sort from low to high, then store eigens and phase velocities.
        """
        eig_val, eig_vec = np.linalg.eigh(self.christoffel)
        args = np.argsort(eig_val)
        eig_val = eig_val[args]
        eig_vec = eig_vec.T[args]

        self._eigenval = eig_val
        self._eigenvec = eig_vec
        self._phase_velocity = np.sign(eig_val)*np.sqrt(np.absolute(eig_val))

    def set_grad_mat(self):
        """
        Calculate the gradient of the Christoffel matrix.
        d/dx_n M_ij = sum_k q_k * ( C_inkj + C_iknj )
        gradmat[n][i][j] =  d/dx_n M_ij (note the indices)
        """
        q = self.direction
        tens = self.stiffness
        gradmat = np.dot(q, tens + np.transpose(tens, (0, 2, 1, 3)))
        gradmat = np.transpose(gradmat, (1, 0, 2))
        self._grad_mat = gradmat

    def set_group_velocity(self):
        """
        Calculate group velocities as the gradient of the phase velocities.
        Powerflow angles are also calculated and stored.
        """
        phase_vel = self.get_phase_velocity()
        eig_vec = self.get_eigenvec()
        gradmat = self.get_grad_mat()

        grad_eig = np.empty((3, 3))
        group_vel = np.empty((3, 3))
        self._group_abs = np.empty(3)
        self._group_dir = np.empty((3, 3))
        self.group_theta = np.empty(3)
        self.group_phi = np.empty(3)
        for pol in xrange(3):
            for cart in xrange(3):
                grad_eig[pol][cart] = \
                np.dot(eig_vec[pol], np.dot(gradmat[cart], eig_vec[pol]))
                # Eigenvalues are the square of the velocity
                # dv/dq = dv^2/dq / (2v)
                group_vel[pol][cart] = grad_eig[pol][cart] / (2*phase_vel[pol])
            self._group_abs[pol] = norm(group_vel[pol])
            self._group_dir[pol] = group_vel[pol] / self._group_abs[pol]

            x = self._group_dir[pol][0]
            z = self._group_dir[pol][2]
            if z >= 1.0-1e-10 or z <= -1.0+1e-10:
                self.group_theta[pol] = 0.0
                self.group_phi[pol] = 0.0
            else:
                self.group_theta[pol] = np.arccos(z)
                sin_theta = np.sqrt(1 - z**2)
                if abs(x) > sin_theta:
                    self.group_phi[pol] = (1.0 - np.sign(x))*0.5*np.pi
                else:
                    self.group_phi[pol] = np.arccos(x/sin_theta)
                if self._group_dir[pol][1] < 0.0:
                    self.group_phi[pol] = 2*np.pi - self.group_phi[pol]
        # In case things go wrong, check if phase_vel == np.dot(group_vel, q)
        self._grad_eig_val = grad_eig
        self._group_velocity = group_vel
        self._cos_pf_angle = np.dot(self._group_dir, self.direction)
        self._powflow_angle = np.arccos(np.around(self._cos_pf_angle, 10))

    def set_hessian_eig(self):
        """
        Calculate the hessian of the eigenvalues.
        Hessian[n][i][j] = d^2 lambda_n / dx_i dx_j
        """
        dynmat = self.christoffel
        eig_val = self.get_eigenval()
        eig_vec = self.get_eigenvec()
        gradmat = self.get_grad_mat()
        hess_mat = self.get_hessian_mat()
       
        hessian = np.zeros((3, 3, 3))
        for n in xrange(3):
            hessian[n] += np.dot(np.dot(hess_mat, eig_vec[n]), eig_vec[n])
            pseudoinv = np.linalg.pinv(eig_val[n]*idmat - dynmat, rcond=1e-10)
            deriv_vec = np.dot(gradmat, eig_vec[n])
            hessian[n] += 2.0 * np.dot(np.dot(deriv_vec, pseudoinv), deriv_vec.T)
            #Take deriv of eigenvec into account: 2 * (d/dx s_i) * pinv_ij * (d_dy s_j)
        self._hessian_eig = hessian
        
    def set_enhancement(self):
        """
        Determine the enhancement factors.
        """
        hessian = self.get_hessian_eig()
        phase_vel = self.get_phase_velocity()
        group_vel = self.get_group_velocity()
        group_abs = self.get_group_abs()

        grad_group = np.empty((3, 3, 3))
        enhance = np.empty(3)
        for n in xrange(3):
            grad_group[n] = hessian[n] / group_abs[n]
            grad_group[n] -= np.outer(group_vel[n], np.dot(hessian[n], group_vel[n])) / (group_abs[n]**3)
            grad_group[n] /= 2.0*phase_vel[n] #grad lambda = 2 * v_p * v_g

            enhance[n] = 1.0 / norm(np.dot(cofactor(grad_group[n]), self.direction))
        self._enhancement = enhance

    def set_enhancement_approx(self, num_steps=8, delta=1e-5):
        """
        Determine the enhancement factors according to a numerical scheme.
        The surface areas of a set of triangles in phase and group space are
        calculated and divided. This is significantly slower and less accurate
        than the analytical approach, but will provide a physically relevant
        value when the enhancement factor is ill defined.

        The surface area is a polygon of n sides where n is num_steps.
        The radius of this polygon is determined by delta, which determines the
        change in theta and phi coordinates relative to the central position.
        """
        phase_grid = np.empty((num_steps+1, 3))
        group_grid = np.empty((num_steps+1, 3, 3))

        center_theta = self.theta
        center_phi = self.phi
        phase_center = self.direction

        for i in xrange(num_steps):
            angle = i*2.0*np.pi/num_steps
            self.set_direction_spherical(center_theta + np.sin(angle)*delta, center_phi + np.cos(angle)*delta)
            phase_grid[i] = self.direction
            group_grid[i] = self.get_group_dir()

        phase_grid[num_steps] = phase_grid[0]
        group_grid[num_steps] = group_grid[0]

        self.set_direction_cartesian(phase_center)
        group_center = self.get_group_dir()

        phase_area = 0.0
        group_area = np.zeros(3)
        tot_angle = np.zeros(3)
        for i in xrange(num_steps):
            phase_area += norm(np.cross(phase_grid[i] - phase_center, phase_grid[i+1] - phase_center))
            for n in xrange(3):
                group_area[n] += norm(np.cross(group_grid[i][n] - group_center[n], group_grid[i+1][n] - group_center[n]))
        self._enhancement = phase_area/group_area

    def find_nopowerflow(self, step_size=0.9, eig_id=2, max_iter=900):
        """
        Attempts to find the closest direction of extremal phase velocity,
        where group and phase directions align. A positive step_size should
        search for maxima, while negative step_size searches for minima.

        Due to the complicated nature of the ray surfaces of the quasi-shear
        modes (eig_id 0 and 1), there is no guarantee that this algorithm
        will converge or reliably find an extremal velocity.

        If a direction has been set already, the search will start from there
        and follow the general direction of power flow. Otherwise, the search
        will start from a randomly chosen point.
        """
        if self.direction is None:
            self.set_direction_random()

        phase_dir = self.direction
        group_dir = self.get_group_dir()

        step_dir = group_dir[eig_id] - phase_dir
        if max_iter <= 0 or norm(step_dir) < 1e-10:
            return
        else:
            self.set_direction_cartesian(phase_dir + step_size*step_dir)
            max_iter -= 1
            self.find_nopowerflow(step_size, eig_id, max_iter)


def de_voigt(C_ij):
    """Turn a 6x6 matrix into a 3x3x3x3 tensor according to Voigt notation."""
    C_ijkl = [[[[C_ij[VOIGT[10*i+j]][VOIGT[10*k+l]]
                 for i in range(3)] for j in range(3)]
                 for k in range(3)] for l in range(3)]
    return C_ijkl

def de_voigt2(vec):
    """Turn a 6-dim vector into a 3x3 tensor according to Voigt notation."""
    T_ij = [[vec[VOIGT[10*i+j]] for i in range(3)] for j in range(3)]
    return T_ij

def hessian_christoffelmat(C):
    """
    Return the hessian of the dynamical matrix.
    Due to the definition of the dynmat (q.C.q), this is independent of q.
    hessianmat[i][j][k][l] = d^2 M_kl / dx_i dx_j (note the indices).
    """
    hessianmat = np.empty((3, 3, 3, 3))
    for i in xrange(3):
        for j in xrange(3):
            for k in xrange(3):
                for l in xrange(3):
                    hessianmat[i][j][k][l] = C[k][i][j][l] + C[k][j][i][l]
    return hessianmat

def get_bulk(C):
    """Return Bulk modulus from stiffness matrix C."""
    return (C[0][0]+C[1][1]+C[2][2] + 2*(C[0][1]+C[0][2]+C[1][2]))/9
def get_shear(C):
    """Return Shear modulus from stiffness matrix C."""
    return ((C[0][0]+C[1][1]+C[2][2]) - (C[0][1]+C[0][2]+C[1][2]) \
            + 3*(C[3][3]+C[4][4]+C[5][5]))/15

def isotropic_velocities(bulk, shear, dens):
    """
    Return primary and secondary sound velocities for an isotropic material.
    Bulk and Shear modulus are assumed to be in GPa, the density in kg/m^3.
    The velocities are returned in km/s.
    """
    primary = np.sqrt(1000.0*(bulk + 4.0*shear/3)/dens)
    secondary = np.sqrt(1000.0*shear/dens)
    return primary, secondary

def get_rot_mat(vector1, vector2):
    """Return a rotation matrix that rotates vector2 towards vector1."""
    vector1 = np.array(vector1)/norm(vector1)
    vector2 = np.array(vector2)/norm(vector2)
    rotvec = np.cross(vector2, vector1)

    sin_angle = norm(rotvec)
    cos_angle = np.sqrt(1.0 - sin_angle*sin_angle)
    if sin_angle > 1e-10:
        dir_vec = rotvec/sin_angle
    else:
        return idmat

    ddt = np.outer(dir_vec, dir_vec)
    skew = np.array([[        0.0, -dir_vec[2],  dir_vec[1]],
                     [ dir_vec[2],         0.0, -dir_vec[0]],
                     [-dir_vec[1],  dir_vec[0],        0.0]])

    mtx = ddt + cos_angle * (idmat - ddt) - sin_angle * skew
    return mtx

def invert_file(filename, theta_column=None, phi_column=None, cart_columns=[]):
    """
    Since the Christoffel tensor is symmetric under inversion, it is only
    necessary to produce half of all data, regardless of crystal symmetry.
    This function will double the data according to inversion symmetry.
    Data will be duplicated from the bottom to the top, and blank lines
    will be reproduced as well.
    
    Keyword arguments:
    filename -- File which contains half of the data.
    theta_column -- Column containing the polar angle.
        Theta -> pi - Theta
    phi_column -- Column containing the azimuthal angle.
        Phi -> Phi +/- pi (in [0, 2pi[)
    cart_columns -- List of columns containing data in cartesian coordinates.
        X -> -X
    """
    infile = open(filename, 'r')
    data = infile.readlines()
    infile.close()

    outfile = open(filename, 'a')
    for linenumber in xrange(len(data)-1, -1, -1):
        line = data[linenumber]
        if line[0] == '#':
            continue
        if line[0] == '\n':
            outfile.write('\n')
            continue
        line = line.split()

        if theta_column is not None:
            line[theta_column] = str( np.pi - float(line[theta_column]) )
        if phi_column is not None:
            phi = float(line[phi_column])
            if phi > np.pi:
                phi = phi - np.pi
            else:
                phi = phi + np.pi
            line[phi_column] = str(phi)
        for column in cart_columns:
            line[column] = str(-float(line[column]))
        line = '\t'.join(line) + '\n'
        outfile.write(line)
    outfile.close()

def determinant(m):
    """Return the determinant of a 3x3 matrix."""
    return (m[0][0] * m[1][1] * m[2][2] -
            m[0][0] * m[1][2] * m[2][1] +
            m[0][1] * m[1][2] * m[2][0] -
            m[0][1] * m[1][0] * m[2][2] +
            m[0][2] * m[1][0] * m[2][1] -
            m[0][2] * m[1][1] * m[2][0])

def norm(v):
    """Return the Pythagorean norm of a 3-dim vector."""
    return np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def cofactor(m):
    """
    Return the cofactor matrix of a 3x3 matrix.
    """
    cof = np.empty((3, 3))

    cof[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1]
    cof[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2]
    cof[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0]

    cof[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2]
    cof[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0]
    cof[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1]

    cof[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1]
    cof[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2]
    cof[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0]

    return cof

