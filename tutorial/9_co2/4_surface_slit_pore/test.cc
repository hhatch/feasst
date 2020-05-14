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

/*
Note:
- the number of Ewald lattice vectors are the same in x,y,z (bad!)
 */

// Define a new Analyze class for error checks
class AnalyzeConfinedFluid_TC34Slit : public feasst::Analyze {
 public:
  AnalyzeConfinedFluid_TC34Slit(shared_ptr<feasst::Pair> pair) : Analyze(pair) {}
  feasst::Accumulator nFluid, zVal;
  void update(const int iMacro) {
    nFluid.accumulate(iMacro);
    for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
      zVal.accumulate(space()->xMol(iMol, 2));
    }
  }
};

int main() {  // CO2, CONFINEMENT_SLIT_PORE
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
                                upper - pairCO2->sig(itype)/2.,
                                lower + pairCO2->sig(itype)/2.,
                                confineDim);
  }
  feasst::PairBarriers pairWall(space.get(), &specbarrier);

  // Create PairHybrid object which encompasses both LJ and wall interaction
  auto pair = feasst::makePairHybrid(space.get());
  pair->addPair(pairCO2.get());
  pair->addPair(&pairWall);
  pair->initEnergy();

  // Initialize acceptance criteria and flat-histogram
  const double activ = exp(-11);
  const double temp = 450; // kelvin
  const int nMolMax = 200;
  auto criteria = feasst::makeCriteriaWLTMMC(
   {{"beta", feasst::str(1./(temp*feasst::idealGasConstant/1e3))},
    {"activ", feasst::str(activ)},
    {"mType", "nmol"},
    {"nMax", feasst::str(nMolMax)}});
  criteria->collectInit();  // initialize updating of collection matrix
  criteria->tmmcInit();     // initialize use of collection matrx (e.g., TMMC)

  // Initialize WLTMMC and trial moves
  feasst::WLTMMC mc(space, pair, criteria);

  mc.weight = 0.4;
  feasst::transformTrial(&mc, "translate");
  feasst::transformTrial(&mc, "rotate");
  
  auto tdel = make_shared<feasst::TrialDelete>(space->addMolListType(0).c_str());
  // tdel->numFirstBeads(10);
  mc.weight = 0.1;
  mc.initTrial(tdel);

  auto tadd = make_shared<feasst::TrialAdd>(space->addMolListType(0).c_str());
  // tadd->numFirstBeads(10);
  mc.weight = 0.1;
  mc.initTrial(tadd);

  // Note that confine can't be called until the space object pointer is
  // initialized. This means it must follow "initTrial" in this example.
  tdel->confine(upper, lower, confineDim);
  tadd->confine(upper, lower, confineDim);

  // Initialize peridoic outputs and checks
  const int nPrint = 1e4;
  mc.initLog("log", nPrint);
  mc.initMovie("movie", nPrint);
  // mc.initXTC("traj", nPrint);
  mc.initColMat("colMat", nPrint);
  mc.initRestart("tmp/rst", nPrint);
  mc.setNFreqCheckE(10*nPrint, 1e-4);
  mc.setNFreqTune(nPrint);

  // Initialize the custom analysis
  auto an = make_shared<AnalyzeConfinedFluid_TC34Slit>(pairCO2);
  an->initFreq(1);    // frequency that Analyze::update() is called
  mc.initAnalyze(an);

  // Initialize parallel windows by dividing order parameter, nmol
  //mc.initWindows(1);
  //mc.runNumSweeps(10, -1);

  // Serial run
  mc.initProduction();
  mc.runNumTrials(1e7);

  ASSERT((an->nFluid.average() > 3) && (an->nFluid.average() < nMolMax - 5),
    "Average number of fluid particles: " << an->nFluid.average());
  ASSERT(an->zVal.min() >= lower, "zVal minimum: " << an->zVal.min());
  ASSERT(an->zVal.max() >= lower, "zVal maximum: " << an->zVal.max());
}

