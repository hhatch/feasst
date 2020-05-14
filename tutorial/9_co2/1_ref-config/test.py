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

class TestCO2_TRAPPE_REFCONF(unittest.TestCase):
  def test(self):
    space = feasst.makeSpace(feasst.args(
      {"dimen" : "3",
       "boxLength" : "43.6254267090798"}))
    pair = feasst.makePairLJCoulEwald(space, feasst.args(
      {"rCut" : "15.",
       "molTypeInForcefield" : "data.CO2_trappe",
       "alphaL" : "5.6",
       "k2max" : "27"}))
    pair.initAtomCut(1)
    pair.equateRcutForAllTypes()

    # read in configuration
    ss = space.install_dir() + "/tutorial/9_co2/1_ref-config/config.xyz"
    space.readXYZAlt(ss)
    pair.addPart()
    pair.initEnergy()

    self.assertAlmostEqual(pair.peTot(), -6742.8401140808373, 15)
    self.assertAlmostEqual(pair.peLJ(), -5351.749077445852, 15)
    self.assertAlmostEqual(pair.peLRC(), -95.224265138825444, 15)
    self.assertAlmostEqual(pair.peQReal(), -1289.8802189291141, 15)
    self.assertAlmostEqual(pair.peQFrr(), 8.1795923734787532, 15)
    self.assertAlmostEqual(pair.peQFrrSelf(), 14.166144940524578, 15)

if __name__ == "__main__":
    unittest.main()
