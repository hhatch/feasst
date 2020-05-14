"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http:#pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

import math
import feasst

nMolMax = 24
feasst.ranInitByDate()
space = feasst.makeSpace(feasst.args(
 {"dimen": "2",
  "boxLength": "40."}))
pair = feasst.makePairPatchKFMulti(space, feasst.args(
  {"rCut": "1.05",  # Note: rCut is the patch well center-center distance
   "patchAngle": "65"}))
# Note: data file sets the bond length, epsilon (depth) and sigma
# Pair Coeffs are in order [type] [eps] [sig]
# Patch unit vectors have 0 sigma
# Center of mass rotation centers of 0 eps (only in cg2_2patch)
ss = space.install_dir() + "/forcefield/data.cg1_2patch_hetero_linear"
pair.initData(ss)

# If you change the sigmas, rCut is not set dependent on sigma, so this
# initIJ may need to be changed if you want rcut-sigma=constant
pair.initIJ()

criteria = feasst.makeCriteriaWLTMMC(feasst.args(
  {"beta": str(1./0.1),   # beta=1/kT
   "activ": str(math.exp(-7)),  # lnz = beta*mu -(debroglie)
   "mType": "nmol",
   "nMax": str(nMolMax),
   "nMin": "0"
  }))
criteria.collectInit()
criteria.tmmcInit()
mc = feasst.WLTMMC(pair, criteria)
mc.weight = 0.4
feasst.addTrialTransform(mc, feasst.args(
  {"transType": "translate",
   "maxMoveParam": "0.1"}))
feasst.addTrialTransform(mc, feasst.args(
  {"transType": "rotate",
   "maxMoveParam": "0.1"}))
mc.weight = 0.2
feasst.insertDeleteTrial(mc, ss)

# begin cluster moves
mc.weight = 1./float(nMolMax)
feasst.gcaTrial(mc)
space.addTypeForCluster(0)

# Set the size of "premicellar aggregates"
# This will print the concentration of "clusters" which are of size <= this value
space.preMicellarAgg(5)

mc.weight = 1./5./float(nMolMax)
feasst.clusterTrial(mc, "clustertrans")
feasst.clusterTrial(mc, "clusterrotate")

# Note that the movie file outputs "patches" as inscribed spheres.
# From PairPatchKFMulti::printxyz line 179:
# make a patch by inscribing sphere of radius r2 inside bead, radius r1=sig
#  distance between center of sphere and bead is a
#  for given patch angle, use law of cosines to derive r2
#  constraint r1 + eps = r2 + a, where eps is small, to avoid clipping
nPrint = int(1e4)
mc.initLog("log", nPrint)
mc.initColMat("colMat", nPrint)
mc.initMovie("movie", nPrint)
mc.initRestart("tmp/rst", nPrint)
mc.setNFreqCheckE(nPrint, 1e-6)
mc.setNFreqTune(nPrint)

# reset cluster statistics after equilibration
space.clusterReset()

mc.runNumSweeps(10, int(1e6))

