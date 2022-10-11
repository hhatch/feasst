# Tutorial in development.
# Compute the interactions between two rigid domains of a protein.

import numpy as np
from pyfeasst import coarse_grain_pdb
from scipy.spatial.transform import Rotation as R

pdb_file = "../../../pyfeasst/tests/1igt.pdb"
chains = {'fc': {'B': range(248, 475), 'D': range(248, 475)}}
fc = coarse_grain_pdb.subset(pdb_file=pdb_file, chains=chains['fc'])
r_com_fc = coarse_grain_pdb.center_of_mass(fc)
print(fc['x_coord'] - r_com_fc[0])

# test that the Euler convention is the same as the one used in FEASST.
#x = [0.6, 0.6, 0.6]
#print(R.from_euler('ZXZ', [90, 0, 0], degrees=True).apply(x))
#print(R.from_euler('ZXZ', [0, 90, 0], degrees=True).apply(x))
#print(R.from_euler('ZXZ', [0, 0, 90], degrees=True).apply(x))
#print(R.from_euler('ZXZ', [180, 0, 0], degrees=True).apply(x))
#print(R.from_euler('ZXZ', [30, 30, 30], degrees=True).apply(x))


