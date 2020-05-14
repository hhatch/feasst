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

class TestCO2_REFCONFIG_CO2_MFI_SILICATE(unittest.TestCase):
  def test(self):

    print("******************************* NOTE *******************************")
    print("* Beware the subtle difference in xyz file formats.")
    print("* FEASST uses the second line as [order_param, lx, ly, lz]")
    print("* Dan uses the second line as [lx, ly, lz]")
    print("* MFI_replicate.xyz follows Dans format,")
    print("* while MFI_replicate_out.xyz from process.py follows FEASSTs format")
    print("* co2.xyz and the input/outputs of concatenate.py follow FEASSTs format")
    print("******************************* NOTE *******************************")

    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "3"}))
    pair = feasst.makePairLJCoulEwald(space, feasst.args(
     {"rCut" : "15",
      "molType" : "none"}))

    # initialize adsorbate, which must be first to use order parameter nMol0
    CO2_data_file = space.install_dir() + "/forcefield/data.CO2_trappe"
    pair.initData(CO2_data_file)

    # add CO2 as a molecule type
    num_CO2 = 20
    for i in range(num_CO2):
      pair.addMol(CO2_data_file)

    # open the json file created by process.py
    with open(space.install_dir() + '/tutorial/9_co2/1_ref-config/MFI_replicate_out_types.json', 'r') as json_file:
      import json
      types = json.load(json_file)

    # initialize framework atom types
    # each atom type should have a corresponding LAMMPS-like data file
    for atom in types:
      data_file_name = space.install_dir() + "/forcefield/data." + atom
      pair.initData(data_file_name)
      for i in range(types[atom]):
        pair.addMol(data_file_name)

    # read the framework + CO2 coordinates from process.py and concatenate.py
    space.readXYZAlt("both.xyz")
    # pair.printXYZ("movie", 1)

    # check if any CO2 has bonds crossing the boundaries
    space.xMolGen()
    for mol in range(0, num_CO2):
      for ox in [1,2]:
        r2 = 0
        ipart = space.mol2part()[mol]
        for dim in range(space.dimen()):
          r2 += (space.x(ipart + ox, dim) - space.x(ipart, dim))**2
        # print(math.sqrt(r2))
        assert(abs(math.sqrt(r2) - 1.16) < 1e-8)

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
    print("lx ", space.boxLength(0))
    print("ly ", space.boxLength(1))
    print("lz ", space.boxLength(2))

    # pair.printXYZ("hi2", 1)

#    mpart = [0]
#    print(pair.multiPartEner(pyfeasst.list_to_int_vector(mpart), 0))
#    print("LJ MFI-CO2", pair.peLJone())
#    print("QReal MFI-CO2", pair.peQRealone())
#    print("nMol", space.nMol())

if __name__ == "__main__":
    unittest.main()
