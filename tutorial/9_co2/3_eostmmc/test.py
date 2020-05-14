"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http:#pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import math
import feasst

openMP = 0
nMolMax = 5
nMolMin = 0
nfreq = int(1e5)
ncfreq = int(1e5)
rCut = 15.
temp = 525
lnz = -8.14
boxl = 30.
molType = "data.CO2_trappe"

# initialize simulation domain
feasst.ranInitByDate()
space = feasst.makeSpace(feasst.args(
  {"dimen" : "3",
   "boxLength" : str(boxl)}))

# initialize pair-wise interactions
pair = feasst.makePairLJCoulEwald(space, feasst.args(
  {"rCut" : str(rCut),
   "molTypeInForcefield" : molType,
   "alphaL" : "5.6",
   "k2max" : "38"}))
pair.initAtomCut(1)
pair.equateRcutForAllTypes()

# acceptance criteria
criteria = feasst.makeCriteriaWLTMMC(feasst.args(
  {"beta" : str(1./(temp*feasst.idealGasConstant/int(1e3))),
   "activ" : str(math.exp(lnz)),
   "mType" : "nmol",
   "nMin" : str(nMolMin),
   "nMax" : str(nMolMax)}))
criteria.collectInit()
criteria.tmmcInit()

# initialize MC simulation object
mc = feasst.WLTMMC(pair, criteria)
mc.weight = 0.4
feasst.transformTrial(mc, "translate")
feasst.transformTrial(mc, "rotate")
mc.weight = 0.1
feasst.insertDeleteTrial(mc, space.addMolListType(0))

# output log, lnpi and movie
mc.initLog("log", nfreq)
mc.initColMat("colMat", ncfreq)
mc.setNFreqCheckE(ncfreq,
                  2e-7)  # absolute (no-percentage) energy tolerance
mc.setNFreqTune(nfreq)
mc.initMovie("movie", nfreq)
# mc.initXTC("xtc", nfreq)
mc.initRestart("tmp/rst", ncfreq)

#production tmmc simulation
if openMP != 0:
  mc.initWindows(1.75,  # exponent that determines size of windows
                 0)    # extra macrostate overlap between processors

mc.runNumSweeps(2,   # number of "sweeps"
               -1)   # maximum number of trials. Infinite if "-1".



