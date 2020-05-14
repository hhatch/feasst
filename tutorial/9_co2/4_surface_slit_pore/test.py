"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http:#pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import unittest
import feasst

"""
Note:
- the number of Ewald lattice vectors are the same in x,y,z (bad!)
"""

class TestCO2_CONFINEMENT_EOSTMMC(unittest.TestCase):
  def test(self):
    feasst.ranInitByDate()
    space = feasst.makeSpace(feasst.args(
     {"dimen" : "3",
      "boxLength" : "30"}))
    space.initBoxLength(space.boxLength(1) + 60, 2)     # stretch z-dimension PBC box length for slab

    # Initialize LJ interactions
    pairCO2 = feasst.makePairLJCoulEwald(space, feasst.args(
     {"rCut" : "15.",
      "molTypeInForcefield" : "data.CO2_trappe",
      "alphaL" : "5.6",
      "k2max" : "38"}))
    pairCO2.initAtomCut(1)
    pairCO2.equateRcutForAllTypes()
    pairCO2.initEnergy()

    # Create an invisible hard slit pore in the z plane.
    specbarrier = feasst.SpeciesBarriers(space.nParticleTypes())
    upper = space.boxLength(0)/2.
    lower = -upper
    confineDim = 2  # z-dimension
    for itype in range(space.nParticleTypes()):
        specbarrier.addHardSlitPore(itype,
                                    upper - pairCO2.sig(itype)/2.,
                                    lower + pairCO2.sig(itype)/2.,
                                    confineDim);
    pairWall = feasst.PairBarriers(space, specbarrier)

    # Create PairHybrid object which encompasses both LJ and wall interaction
    pair = feasst.makePairHybrid(space)
    pair.addPair(pairCO2)
    pair.addPair(pairWall)
    pair.initEnergy()

    # Initialize acceptance criteria and flat-histogram
    import math
    activ = math.exp(-11)
    temp = 450 # kelvin
    nMolMax = 200
    criteria = feasst.makeCriteriaWLTMMC(feasst.args(
     {"beta" : str(1./(temp*feasst.idealGasConstant/int(1e3))),
      "activ" : str(activ),
      "mType" : "nmol",
      "nMax" : str(nMolMax)}))
    criteria.collectInit()  # initialize updating of collection matrix
    criteria.tmmcInit()     # initialize use of collection matrx (e.g., TMMC)

    # Initialize WLTMMC and trial moves
    mc = feasst.WLTMMC(space, pair, criteria)

    mc.weight = 0.4
    feasst.transformTrial(mc, "translate")
    feasst.transformTrial(mc, "rotate")

    tdel = feasst.TrialDelete(space.addMolListType(0))
    # tdel.numFirstBeads(10)
    mc.weight = 0.1
    mc.initTrial(tdel)

    tadd = feasst.TrialAdd(space.addMolListType(0))
    # tadd.numFirstBeads(10)
    mc.weight = 0.1
    mc.initTrial(tadd)

    # Note that confine can't be called until the space object pointer is
    # initialized. This means it must follow "initTrial" in this example.
    tdel.confine(upper, lower, confineDim)
    tadd.confine(upper, lower, confineDim)

    # Initialize peridoic outputs and checks
    nPrint = int(1e4)
    mc.initLog("log", nPrint)
    mc.initMovie("movie", nPrint)
    # mc.initXTC("traj", nPrint)
    mc.initColMat("colMat", nPrint)
    mc.initRestart("tmp/rst", nPrint)
    mc.setNFreqCheckE(10*nPrint, 1e-4)
    mc.setNFreqTune(nPrint)

    # Initialize parallel windows by dividing order parameter, nmol
    #mc.initWindows(1)
    #mc.runNumSweeps(10, -1)

    # Serial run
    mc.initProduction()
    mc.runNumTrials(int(1e7))

if __name__ == "__main__":
    unittest.main()
