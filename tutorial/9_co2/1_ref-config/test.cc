/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "feasst.h"

void expect_near91(const double val1, const double val2, const double tolerance, const char* comment) {
  const double diff = val1 - val2;
  if (fabs(diff) > tolerance) {
    ASSERT(0, "Expected " << MAX_PRECISION << val1 << " to be the same as " << val2 << " within tolerance " << tolerance << " for case " << comment);
  }
}

int main() {  // CO2_TRAPPE, REFCONF
  auto space = feasst::makeSpace(
    {{"dimen", "3"},
     {"boxLength", "43.6254267090798"}});
  auto pair = feasst::makePairLJCoulEwald(space,
    {{"rCut", "15."},
     {"molTypeInForcefield", "data.CO2_trappe"},
     {"alphaL", "5.6"},
     {"k2max", "27"}});
  pair->initAtomCut(1);
  pair->equateRcutForAllTypes();

  // read in configuration
  std::stringstream ss;
  ss << space->install_dir() << "/tutorial/9_co2/1_ref-config/config.xyz";
  std::ifstream file(ss.str().c_str());
  pair->readXYZ(file);
  pair->initEnergy();

  expect_near91(pair->peTot(), -6742.8401140808373, feasst::DTOL, "Total");
  expect_near91(pair->peLJ(), -5351.749077445852, feasst::DTOL, "LJ");
  expect_near91(pair->peLRC(), -95.224265138825444, feasst::DTOL, "LRC");
  expect_near91(pair->peQReal(), -1289.8802189291141, feasst::DTOL, "QReal");
  expect_near91(pair->peQFrr(), 8.1795923734787532, feasst::DTOL, "QFrr");
  expect_near91(pair->peQFrrSelf(), 14.166144940524578, feasst::DTOL, "QFrrSelf");
}
