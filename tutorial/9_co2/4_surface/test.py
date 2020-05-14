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

def process(in_file_name):
  out_file_name = in_file_name + "_out"

  with open(in_file_name + ".xyz", "r") as f:
    content = f.readlines()

  out_file = open(out_file_name + ".xyz", "w")
  print(content[0], file=out_file, end='')
  print("-1 " + content[1], file=out_file, end='')

  data = {}
  atomTypes = list()
  # numTypes = list()
  for line in range(len(content)):
    if line >= 2:
      name = content[line].split(" ")[0]
      # print("name", name)
      if name not in atomTypes:
        atomTypes.append(name)
        data[name] = 1
      else:
        data[name] = data[name] + 1

  for atom in atomTypes:
    for line in range(len(content)):
      if line >= 2:
        name = content[line].split(" ")[0]
        if name == atom:
          print(content[line], file=out_file, end='')

  import json
  with open(out_file_name + "_types.json", 'w') as json_file:
    json.dump(data, json_file)

  out_file.close()

class TestCO2_MFI_SILICATE(unittest.TestCase):
  def test(self):

    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "3"}))
    in_file_name = space.install_dir() + '/tutorial/9_co2/4_surface/MFI_replicate'
    process(in_file_name)
    pair = feasst.makePairLJCoulEwald(space, feasst.args(
     {"rCut" : "10",
      "molType" : "none"}))

    # open the json file created by process.py
    with open(in_file_name + "_out_types.json", 'r') as json_file:
      import json
      types = json.load(json_file)

    # initialize adsorbate, which must be first to use order parameter nMol0
    CO2_data_file = space.install_dir() + "/forcefield/data.CO2_trappe"
    pair.initData(CO2_data_file)

    # initialize framework atom types
    # each atom type should have a corresponding data. LAMMPS-formatted file
    for atom in types:
      data_file_name = space.install_dir() + "/forcefield/data." + atom
      pair.initData(data_file_name)
      for i in range(types[atom]):
        pair.addMol(data_file_name)

    # read the framework coordinates sorted by process.py
    space.readXYZAlt(in_file_name + "_out.xyz")
    pair.printXYZ("movie", 1)

#    # shift the framework to be centered at the origin
#    xAdd = list()
#    for dim in range(space.dimen()):
#      xAdd.append(-space.boxLength(dim)/2.)
#    for iMol in range(space.nMol()):
#      space.transMol(iMol, pyfeasst.list_to_double_vector(xAdd))

    # initialize LJ cutoff
    #print("should mess up here")
    #pair.cutShift(1);
    #print("should mess up here")
    pair.linearShift(1);

    # initialize Ewald and total energy
    alpha = 5.6
    k2max = 27
    pair.initKSpace(alpha, k2max)
    pair.initEnergy()

    # initialize acceptance criteria
    temp = 500   # kelvin, assumes epsilons are kJ/mol
    nMolMax = 20       # maximum number of CO2 molecules
    activ = math.exp(-6)   # activity,z of CO2
    criteria = feasst.makeCriteriaWLTMMC(feasst.args(
     {"beta" : str(1./(temp*feasst.idealGasConstant/1e3)),
      "activ" : str(activ),
      "mType" : "nmol0", # order parameter is the number of molecules of adsorbate
      "nMax" : str(nMolMax)}))
    for atom in types:
      criteria.addActivity(math.exp(-1))   # activity of framework doesn't matter, but error checked
    criteria.collectInit()
    criteria.tmmcInit()
    mc = feasst.WLTMMC(pair, criteria)

    # initialize trial moves
    ttrans = feasst.TrialTransform("translate")
    ttrans.selectType(CO2_data_file)
    mc.weight = 0.4
    mc.initTrial(ttrans)

    trot = feasst.TrialTransform("rotate")
    trot.selectType(CO2_data_file)
    mc.weight = 0.4
    mc.initTrial(trot)

    mc.weight = 0.1
    feasst.insertDeleteTrial(mc, CO2_data_file)

    # initialize peridoic outputs and checks
    nPrint = int(1e4)
    mc.initLog("log", nPrint)
    mc.initMovie("movie", nPrint)
    # mc.initXTC("traj", nPrint)
    mc.initColMat("colMat", nPrint)
    mc.initRestart("tmp/rst", nPrint)
    mc.setNFreqCheckE(10*nPrint, 1e-4)
    mc.setNFreqTune(nPrint)

    # run serial
    mc.initProduction()
    for sweeps in range(1, 2):
      mc.runNumSweeps(sweeps)   # note: if using OMP, runNumSweeps starts over
      criteria.writeRestart("colMat_sweep" + str(criteria.nSweep()) +
                                  "_trial" + str(mc.nAttempts()))

if __name__ == "__main__":
    unittest.main()
