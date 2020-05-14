"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import math
import unittest
import feasst
import pyfeasst

class TestCO2_REFCONFIG_CO2_CO2(unittest.TestCase):
  def test(self):

    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "3"}))
    pair = feasst.makePairLJCoulEwald(space, feasst.args(
     {"rCut" : "15",
      "molType" : "none"}))

    # initialize adsorbate, which must be first to use order parameter nMol0
    CO2_data_file = space.install_dir() + "/forcefield/data.CO2_trappe"
    pair.initData(CO2_data_file)

    # add 20 CO2
    for i in range(20):
      pair.addMol(CO2_data_file)

    # read the CO2 coordinates
    space.readXYZAlt("co2.xyz")
    pair.printXYZ("movie", 1)

    # loop through the CO2 and see which ones have bonds crossing the boundaries
    space.xMolGen()
    for mol in range(0, space.nMol()):
      for ox in [1,2]:
        r2 = 0
        ipart = space.mol2part()[mol]
        for dim in range(space.dimen()):
          r2 += (space.x(ipart + ox, dim) - space.x(ipart, dim))**2
        # print(math.sqrt(r2))
        assert(abs(math.sqrt(r2) - 1.16) < 1e-8)

#    # shift the framework to be centered at the origin
#    xAdd = list()
#    for dim in range(space.dimen()):
#      xAdd.append(-space.boxLength(dim)/2.)
#    for iMol in range(space.nMol()):
#      space.transMol(iMol, pyfeasst.list_to_double_vector(xAdd))

    # initialize LJ cutoff
    pair.initAtomCut(1)
    pair.equateRcutForAllTypes()
    pair.linearShift(1);

    # initialize Ewald and total energy
    alpha = 5.6
    k2max = 27
    pair.initKSpace(alpha, k2max)
    pair.initEnergy()

    conversion = 1./4.184
    print("Total(kcal/mol):", pair.peTot()*conversion)
    print("LJ(kcal/mol):", pair.peLJ()*conversion)
    print("QReal(kcal/mol):", pair.peQReal()*conversion)
    print("QFrr(kcal/mol):", pair.peQFrr()*conversion)
    print("QFrrSelf(kcal/mol):", pair.peQFrrSelf()*conversion)

    pair.printXYZ("hi2", 1)

    mpart = [0]
    print(pair.multiPartEner(pyfeasst.list_to_int_vector(mpart), 0))
    print("Interaction of the first site/atom with the system")
    print("LJ MFI-CO2", pair.peLJone())
    print("QReal MFI-CO2", pair.peQRealone())
    print("nMol", space.nMol())

if __name__ == "__main__":
    unittest.main()
