/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include "pair_barriers.h"
#include "pair_hybrid.h"
#include "pair_lj.h"

TEST(PairBarriers, one_component_hardwall) {
  feasst::Space space(3);
  space.initBoxLength(8, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);
  std::string addMolType("../forcefield/data.lj");
  space.addMolInit(addMolType);
  for (int iMol = 0; iMol < 50; ++iMol) {
    space.addMol(addMolType);
  }
  space.initBoxLength(30, 0);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addHardSlitPore(0, 5, -5, 0);
  feasst::PairBarriers pair(&space, &specbarrier);

  pair.initEnergy();
  EXPECT_NEAR(0, pair.peTot(), 1e-14);
  pair.printXYZ("tmp/wall5", 1);

  for (int iMol = 0; iMol < 50; ++iMol) {
    space.addMol(addMolType);
  }
  pair.initEnergy();
  EXPECT_LT(NUM_INF, pair.peTot());
}

TEST(PairBarriers, two_component_hardwall) {
  feasst::Space space(3);
  space.initBoxLength(10, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);
  std::string addMolTypeB("../forcefield/data.ljb");
  space.addMolInit(addMolTypeB);

  EXPECT_EQ (space.nParticleTypes(), 2);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addHardSlitPore(0, 5, -5, 0);
  specbarrier.addHardSlitPore(1, 5, -5, 0);

  feasst::PairBarriers pair(&space, &specbarrier);

  for (int iMol = 0; iMol < 25; ++iMol) {
    space.addMol(addMolTypeA);
  space.addMol(addMolTypeB);
  }
  space.initBoxLength(30, 0);

  pair.initEnergy();
  EXPECT_NEAR(0, pair.peTot(), 1e-14);
  pair.printXYZ("tmp/wall6", 1);

  for (int iMol = 0; iMol < 25; ++iMol) {
    space.addMol(addMolTypeA);
  space.addMol(addMolTypeB);
  }

  pair.initEnergy();
  EXPECT_LT(NUM_INF, pair.peTot());
}

TEST(PairBarriers, two_component_diff_hardwall) {
  feasst::Space space(3);
  space.initBoxLength(10, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);
  std::string addMolTypeB("../forcefield/data.ljb");
  space.addMolInit(addMolTypeB);

  EXPECT_EQ (space.nParticleTypes(), 2);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addHardSlitPore(0, 5, -5, 0);
  specbarrier.addHardSlitPore(1, 8, -8, 0);

  feasst::PairBarriers pair(&space, &specbarrier);

  for (int iMol = 0; iMol < 250; ++iMol) {
    space.addMol(addMolTypeA);
  space.addMol(addMolTypeB);
  }
  space.initBoxLength(30, 0);

  pair.initEnergy();
  EXPECT_NEAR(0, pair.peTot(), 1e-14);
  space.printXYZ("tmp/wall7", 1);

  for (int iMol = 0; iMol < 250; ++iMol) {
    space.addMol(addMolTypeA);
  space.addMol(addMolTypeB);
  }

  pair.initEnergy();
  EXPECT_LT(NUM_INF, pair.peTot());
}

TEST(PairBarriers, one_component_sqwwall) {
  feasst::Space space(3);
  space.initBoxLength(10, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);

  EXPECT_EQ (space.nParticleTypes(), 1);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addSqwSlitPore(0, 5, -5, 1.234, 5.0, 0);

  feasst::PairBarriers pair(&space, &specbarrier);

  for (int iMol = 0; iMol < 250; ++iMol) {
    space.addMol(addMolTypeA);
  }
  space.initBoxLength(30, 0);

  pair.initEnergy();
  EXPECT_NEAR(250*-1.234, pair.peTot(), 1.0e-9);
  pair.writeRestart("tmp/barrst");
}

TEST(PairBarriers, two_component_sqwwall) {
  feasst::Space space(3);
  space.initBoxLength(10, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);
  std::string addMolTypeB("../forcefield/data.ljb");
  space.addMolInit(addMolTypeB);

  EXPECT_EQ (space.nParticleTypes(), 2);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addSqwSlitPore(0, 5, -5, 1.234, 5.0, 0);
  specbarrier.addSqwSlitPore(1, 5, -5, 2.345, 5.0, 0);

  feasst::PairBarriers pair(&space, &specbarrier);

  for (int iMol = 0; iMol < 250; ++iMol) {
    space.addMol(addMolTypeA);
    space.addMol(addMolTypeB);

  }
  space.initBoxLength(30, 0);

  pair.initEnergy();
  EXPECT_NEAR(250*-1.234 + 250*-2.345, pair.peTot(), 1.0e-9);
}

TEST(PairBarriers, one_component_sqwcylinderX) {
  feasst::Space space(3);
  space.initBoxLength(30, 0);
  space.initBoxLength(7, 1);
  space.initBoxLength(7, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);

  EXPECT_EQ (space.nParticleTypes(), 1);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  const vector < double > pt (3, 0);
  specbarrier.addSqwCylinder (0, pt, 5.1, 1.234, 5.1, 0);

  feasst::PairBarriers pair(&space, &specbarrier);

  // Box purposefully designed to be too small so particles are created in a confined space
  EXPECT_FALSE (specbarrier.safe(pair));

  // Cylinder "fills" space with interactions
  int npart = 300;
  for (int iMol = 0; iMol < npart; ++iMol) {
    space.addMol(addMolTypeA);
  }
  pair.initEnergy();
  EXPECT_NEAR(npart*-1.234, pair.peTot(), 1.0e-9);

  // Expand the box size to "safe" level
  space.initBoxLength(30, 1);
  space.initBoxLength(30, 2);
  EXPECT_TRUE (specbarrier.safe(pair));
}

TEST(PairBarriers, one_component_sqwcylinderY) {
  feasst::Space space(3);
  space.initBoxLength(7, 0);
  space.initBoxLength(30, 1);
  space.initBoxLength(7, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);

  EXPECT_EQ (space.nParticleTypes(), 1);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  const vector < double > pt (3, 0);
  specbarrier.addSqwCylinder (0, pt, 5.1, 1.234, 5.1, 1);

  feasst::PairBarriers pair(&space, &specbarrier);

  // Box purposefully designed to be too small so particles are created in a confined space
  EXPECT_FALSE (specbarrier.safe(pair));

  // Cylinder "fills" space with interactions
  int npart = 300;
  for (int iMol = 0; iMol < npart; ++iMol) {
    space.addMol(addMolTypeA);
  }
  pair.initEnergy();
  EXPECT_NEAR(npart*-1.234, pair.peTot(), 1.0e-9);

  // Expand the box size to "safe" level
  space.initBoxLength(30, 0);
  space.initBoxLength(30, 2);
  EXPECT_TRUE (specbarrier.safe(pair));
}

TEST(PairBarriers, one_component_sqwcylinderZ) {
  feasst::Space space(3);
  space.initBoxLength(7, 0);
  space.initBoxLength(7, 1);
  space.initBoxLength(30, 2);

  std::string addMolTypeA("../forcefield/data.lj");
  space.addMolInit(addMolTypeA);

  EXPECT_EQ (space.nParticleTypes(), 1);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  const vector < double > pt (3, 0);
  specbarrier.addSqwCylinder (0, pt, 5.1, 1.234, 5.1, 2);

  feasst::PairBarriers pair(&space, &specbarrier);

  // Box purposefully designed to be too small so particles are created in a confined space
  EXPECT_FALSE (specbarrier.safe(pair));

  // Cylinder "fills" space with interactions
  int npart = 300;
  for (int iMol = 0; iMol < npart; ++iMol) {
    space.addMol(addMolTypeA);
  }
  pair.initEnergy();
  EXPECT_NEAR(npart*-1.234, pair.peTot(), 1.0e-9);

  // Expand the box size to "safe" level
  space.initBoxLength(30, 0);
  space.initBoxLength(30, 1);
  EXPECT_TRUE (specbarrier.safe(pair));
}

// define some confinement geometries to pass as arguments
struct pore_args {
  double boxLength, upper0, lower0, hybridCut;
  bool expect_safe;
  std::string type;
};

void testIfSafe(pore_args pore) {
  feasst::Space space(3);
  space.initBoxLength(pore.boxLength);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});
  pairLJ.initData("../forcefield/data.lj");
  pairLJ.initData("../forcefield/data.ljb");
  pairLJ.linearShift(true);

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  if (pore.type == "hardSlit") {
    specbarrier.addHardSlitPore(0, pore.upper0, pore.lower0, 0);
    specbarrier.addHardSlitPore(1, 10, -10, 0);
  } else if (pore.type == "sqwSlit") {
    specbarrier.addSqwSlitPore(0, pore.upper0, pore.lower0, 1.234, 1, 0);
    specbarrier.addSqwSlitPore(1, 10, -10, 1.234, 1, 0);
  } else {
    ASSERT(0, "unrecognized pore.type(" << pore.type << ")");
  }
  feasst::PairBarriers pairBarr(&space, &specbarrier);

  feasst::PairHybrid pair(&space, {{"rCut", feasst::str(pore.hybridCut)}});
  pair.addPair(&pairLJ);
  pair.addPair(&pairBarr);
  pair.initEnergy();

  if (pore.expect_safe) {
    EXPECT_TRUE (specbarrier.safe(pair));
  } else {
    EXPECT_TRUE (!specbarrier.safe(pair));
  }
}

TEST(PairBarriers, isPoreSafe) {
  pore_args pore;
  pore.type.assign("hardSlit");
  pore.boxLength = 30;
  pore.upper0 = 5;
  pore.lower0 = -5;
  pore.hybridCut = 1.234;
  pore.expect_safe = true;
  testIfSafe(pore);

  pore.boxLength = 20;
  pore.expect_safe = false;
  testIfSafe(pore);  // fails for second component

  pore.boxLength = 30;
  pore.upper0 = -10;
  pore.lower0 = -20;
  pore.expect_safe = false;
  testIfSafe(pore);  // fails for second component

  pore.boxLength = 30;
  pore.upper0 = 20;
  pore.lower0 = 10;
  pore.expect_safe = false;
  testIfSafe(pore);  // species 1's wall overlaps boundary

  pore.type.assign("sqwSlit");
  pore.upper0 = 5.;
  pore.lower0 = -5;
  pore.expect_safe = true;
  testIfSafe(pore);

  pore.boxLength = 20.;
  pore.expect_safe = false;
  testIfSafe(pore);

  pore.boxLength = 30.;
  pore.upper0 = -10;
  pore.lower0 = -20;
  pore.expect_safe = false;
  testIfSafe(pore);

  pore.upper0 = 20;
  pore.lower0 = 10;
  pore.expect_safe = false;
  testIfSafe(pore);

  pore.boxLength = 30.;
  pore.upper0 = 10;
  pore.lower0 = -10;
  pore.expect_safe = true;
  pore.hybridCut = 9.999;
  testIfSafe(pore);

  pore.hybridCut = 10.0001; // too big
  pore.expect_safe = false;
  testIfSafe(pore);
}

struct cylPoreArgs {
  vector<double> point;
  bool expect_safe;
  double hybridCut, boxLength, arg;
};

void isCylPoreSafe(cylPoreArgs pore) {
  feasst::Space space(3);
  space.initBoxLength(pore.boxLength);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"cutType", "linearShift"}});

  feasst::SpeciesBarriers specbarrier (space.nParticleTypes());
  specbarrier.addSqwCylinder (0, pore.point, pore.arg, 1.234, pore.arg, 0);
  feasst::PairBarriers pairBarr(&space, &specbarrier);

  feasst::PairHybrid pair(&space, {{"rCut", feasst::str(pore.hybridCut)}});
  pair.addPair(&pairLJ);
  pair.addPair(&pairBarr);

  if (pore.expect_safe) {
    EXPECT_TRUE (specbarrier.safe(pair));
  } else {
    EXPECT_TRUE (!specbarrier.safe(pair));
  }
}

TEST(PairBarriers, sqwCylinder) {
  cylPoreArgs pore;
  pore.point.resize(3);
  pore.expect_safe = true;
  pore.hybridCut = 1.234;
  pore.boxLength = 30;
  pore.arg = 5.1;
  isCylPoreSafe(pore);

  pore.point[2] = 10;
  pore.expect_safe = false;
  isCylPoreSafe(pore);

  pore.point[2] = -10;
  pore.expect_safe = false;
  isCylPoreSafe(pore);

  pore.point[2] = 0;
  pore.hybridCut = 1.9999;
  pore.expect_safe = true;
  isCylPoreSafe(pore);

  pore.boxLength = 20;
  pore.hybridCut = 2.001;
  pore.arg = 9;
  pore.expect_safe = false;
  isCylPoreSafe(pore);
}

TEST(PairBarriers, one_component_ljwall) {
  auto space = feasst::makeSpace({{"dimen", "3"}, {"boxLength", "30"}});
  std::string addMolTypeA("../forcefield/data.lj");
  space->addMolInit(addMolTypeA);
  feasst::SpeciesBarriers specbarrier (space->numSiteTypes());


  specbarrier.addSqwSlitPore(0, 5, -5, 1.0, 1.0, 0);
  specbarrier.specbarriers().back()->features().back()->initLJ(2./15., 9., 3.);
  const double pe_expected = 2*((2./15.)*pow(5, -9) - pow(5, -3));
  feasst::PairBarriers pair(space.get(), &specbarrier);
  std::vector<double> xAdd(space->dimen(), 0.);
  pair.addMol(xAdd, addMolTypeA);
  pair.initEnergy();
  EXPECT_NEAR(pe_expected, pair.peTot(), feasst::DTOL);

  // test with only one LJ surface
  specbarrier.specbarriers(0)->features(0)->initLJ(2./15., 9., 3., 1);
  //specbarrier.specbarriers().back()->features().back()->initLJ(2./15., 9., 3., 1);
  pair.initEnergy();
  EXPECT_NEAR(pe_expected/2., pair.peTot(), feasst::DTOL);
}
