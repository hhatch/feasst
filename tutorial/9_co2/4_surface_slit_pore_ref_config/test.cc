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

int main() {  // CO2, CONFINEMENT_SLIT_PORE_REF_CONFIG
  feasst::ranInitByDate();
  auto space = feasst::makeSpace(
   {{"dimen", "3"},
    {"boxLength", "30"}});
  space->initBoxLength(space->boxLength(1)+60, 2);     // stretch z-dimension PBC box length for slab

  // Initialize LJ interactions
  auto pairCO2 = feasst::makePairLJCoulEwald(space,
   {{"rCut", "15."},
    {"molTypeInForcefield", "data.CO2_trappe"},
    {"alphaL", "5.6"},
    {"k2max", "38"}});
  pairCO2->initAtomCut(1);
  pairCO2->equateRcutForAllTypes();
  pairCO2->initEnergy();

  // Create an invisible hard slit pore in the z plane.
  feasst::SpeciesBarriers specbarrier(space->nParticleTypes());
  const double upper = space->boxLength(0)/2., lower = -upper;
  const int confineDim = 2;  // z-dimension
  for (int itype = 0; itype < space->nParticleTypes(); ++itype) {
    specbarrier.addHardSlitPore(itype,
                                upper - pairCO2->sig(itype)/2,
                                lower + pairCO2->sig(itype)/2,
                                confineDim);
  }
  feasst::PairBarriers pairWall(space.get(), &specbarrier);

  // Create PairHybrid object which encompasses both LJ and wall interaction
  auto pair = feasst::makePairHybrid(space.get());
  pair->addPair(pairCO2.get());
  pair->addPair(&pairWall);
  pair->initEnergy();

  // Add molecule in the center
  std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd);
  pair->addPart();
  pair->initEnergy();
  // pair->printXYZ("hi", 1);
  ASSERT(pair->peTot() < 0, "pe: " << pair->peTot());

  // Move molecule just up to the hard barrier
  xAdd[2] = 15. - pairCO2->sig(1)/2. - 10.*feasst::DTOL - space->x(0, 2);
  space->transMol(0, xAdd);
  pair->initEnergy();
  ASSERT(pair->peTot() < 0, "pe: " << pair->peTot());

  // Move molecule just beyond the hard barrier
  xAdd[2] = 15. - pairCO2->sig(1)/2. + 10.*feasst::DTOL - space->x(0, 2);
  space->transMol(0, xAdd);
  pair->initEnergy();
  // pair->printXYZ("hi", 0);
  ASSERT(pair->peTot() > 1e200, "pe: " << pair->peTot());

  // Rotate the molecule to lie along the z-axis and test again
  space->qMolAlt(0, 0, sqrt(2)/2.);
  space->qMolAlt(0, 1, 0.);
  space->qMolAlt(0, 2, sqrt(2)/2.);
  space->qMolAlt(0, 3, 0.);
  space->quat2pos(0);

  // Move molecule just up to the hard barrier
  xAdd[2] = 15. - 1.16 - pairCO2->sig(1)/2. - 10.*feasst::DTOL - space->x(0, 2);
  space->transMol(0, xAdd);
  pair->initEnergy();
  
  // Move molecule just up to the hard barrier
  xAdd[2] = 15. - 1.06 - pairCO2->sig(1)/2. + 10.*feasst::DTOL - space->x(0, 2);
  space->transMol(0, xAdd);
  pair->initEnergy();
  // space->writeRestart("space");
  // pair->printXYZ("hi", 0);
  ASSERT(pair->peTot() > 1e200, "pe: " << pair->peTot());
}

