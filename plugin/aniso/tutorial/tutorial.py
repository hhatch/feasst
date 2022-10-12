"""
Tutorial in development.
Compute the interactions between two rigid domains of a protein.

Compute the energy between the two domains using AMBER.
If explicit water, this requires equilibration and ensemble average.
Instead, for first pass, start with faster model.
Perhaps not even implicit solvent. Model charges with a screened electrostatic (Yukawa) potential.

Does AMBER have a potential of the form U ~ q_i q_j exp(-kappa r_ij) / r_ij?
where q is the charge, r is the separation distance, and kappa is related to the ionic strength.

First, compute and plot the energy as a function distance, r, along the z-axis with fixed orientation.
This will inform the algorithm used to obtain the r_h and r_c.
Let r_h be the minimum r for U(r_h) < U_{max}.
Let r_c be the 'reasonable' truncation distance that approximates the 'zero' interaction, infinite
separation distance.

Once r_h and r_c are determined algorithmically, grid the energy in 'z'
z = [(r - r_h)/(r_c - r_h)]^gamma
where the choice of gamma optimally stretches r to increase resolution at short distances.

Repeat for each of 5 angle parameters to generate a potential table usable for MC simulations.
For an estimated production quality pair interaction table, there could easily be 41 million
unique configurations to precalculate.
Parallelization and optimization (Python API?) could help.
"""

import numpy as np
from pyfeasst import coarse_grain_pdb
from scipy.spatial.transform import Rotation as R

# select a domain of 1igt
pdb_file = "../../../pyfeasst/tests/1igt.pdb"
domain_name = 'fab1'

chains = {'hinge': {'B': range(236, 244), 'D': range(236, 244)},
          'fc': {'B': range(248, 475), 'D': range(248, 475)},
          'fab1': {'A': range(1, 215), 'B': range(1, 230)},
          'fab2': {'C': range(1, 215), 'D': range(1, 230)},
          'fv1': {'A': range(1, 109), 'B': range(1, 113)},
          'fv2': {'C': range(1, 109), 'D': range(1, 113)},
          'ch1_1': {'A': range(109, 215), 'B': range(113, 230)},
          'ch1_2': {'C': range(109, 215), 'D': range(113, 230)},
          'ch2': {'B': range(248, 361), 'D': range(248, 361)},
          'ch3': {'B': range(361, 475), 'D': range(361, 475)}}
domain = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains[domain_name])

# compute the center of mass of the domain
r_com = coarse_grain_pdb.center_of_mass(domain)

def print_xyz(domain, r_com, file_name, sph=[0, 0, 0], euler=[0., 0., 0.]):
    """
    Print the reference coordinates of the domain to xyz format file.

    First, rotate the domain by the given Euler angles (ZXZ convention)
    This is the "so-called" proper, intrinsic, active rotation rotation convention.
    First, rotate by an angle phi:[-pi, pi] about the z-axis
    Second, rotate by an angle theta:[0, pi] about the new x-axis (intrinsic)
    Third, rotate by an angle psi:[-pi, pi] about the new z-axis
    See https://en.wikipedia.org/wiki/Euler_angles

    Second, translate the domain by the given spherical coordinates.
    Spherical coordinates are given in the mathematics format of
    (radial distance, azimuthal angle[0-2pi], polar angle:[0-pi]).
    """
    fl = open(file_name, 'w')
    fl.write(str(len(domain)) + '\n-1\n')
    if sph[0] < 1e-8:
        new_com = [0, 0, 0]
    else:
        new_com = [sph[0]*np.cos(sph[1])*np.sin(sph[2]),
                   sph[0]*np.sin(sph[1])*np.sin(sph[2]),
                   sph[0]*np.cos(sph[2])]
    for index,_ in enumerate(domain['x_coord']):
        x = domain['x_coord'].values[index] - r_com[0]
        y = domain['y_coord'].values[index] - r_com[1]
        z = domain['z_coord'].values[index] - r_com[2]
        xn = R.from_euler('ZXZ', euler).apply([x, y, z])
        fl.write('0 ' + str(round(xn[0] + new_com[0], 6)) + ' '
                      + str(round(xn[1] + new_com[1], 6)) + ' '
                      + str(round(xn[2] + new_com[2], 6)) + '\n')
    fl.close()

print_xyz(domain, r_com, file_name=domain_name+'_0_0_0_0_0_0.xyz')
print_xyz(domain, r_com, file_name=domain_name+'_50_0_0_0_0_0.xyz', sph=[50, 0, 0], euler=[0, 0, 0])

# Compute energy of interaction between the two xyz file structures using AMBER.
# Does AMBER have an easy-to-use Python API for something like this?

