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
#include "barrier.h"

TEST(Barrier, HardSlitPoreX) {
  feasst::Barrier barrier;
  barrier.addHardSlitPore(5, -5, 0);

  vector<double> x(3, 0.);
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[0] = -4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = -5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, HardSlitPoreY) {
  feasst::Barrier barrier;
  barrier.addHardSlitPore(5, -5, 1);

  vector<double> x(3, 0.);
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[1] = -4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = -5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, HardSlitPoreZ) {
  feasst::Barrier barrier;
  barrier.addHardSlitPore(5, -5, 2);

  vector<double> x(3, 0.);
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[2] = -4.99;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = -5.01;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwSlitPoreX) {
  feasst::Barrier barrier;
  double upper = 8.0, lower = 3.5, eps = 1.234, width = 1.0;
  int dimen = 0;
  barrier.addSqwSlitPore (upper, lower, eps, width, dimen);

  vector<double> x(3, 6.0);

  x[0] = 6.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.5001;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 7.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.49999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 8.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[0] = 3.4999;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwSlitPoreY) {
  feasst::Barrier barrier;
  double upper = 8.0, lower = 3.5, eps = 1.234, width = 1.0;
  int dimen = 1;
  barrier.addSqwSlitPore (upper, lower, eps, width, dimen);

  vector<double> x(3, 6.0);

  x[1] = 6.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.5001;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 7.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.49999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 8.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[1] = 3.4999;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwSlitPoreZ) {
  feasst::Barrier barrier;
  double upper = 8.0, lower = 3.5, eps = 1.234, width = 1.0;
  int dimen = 2;
  barrier.addSqwSlitPore (upper, lower, eps, width, dimen);

  vector<double> x(3, 6.0);

  x[2] = 6.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.5001;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 7.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.49999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 8.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));

  x[2] = 3.4999;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderXY) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 0;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[1] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[1] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderXZ) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 0;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[2] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[2] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderYX) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 1;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[0] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[0] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderYZ) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 1;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[2] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[2] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[2] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[2] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderZX) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 2;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[0] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[0] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[0] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[0] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}

TEST(Barrier, SqwCylinderZY) {
  feasst::Barrier barrier;
  double radius = 5.0, eps = 1.234, width = 1.0;
  int dimen = 2;
  vector < double > pt(3, 0);
  barrier.addSqwCylinder (pt, radius, eps, width, dimen);

  vector<double> x(3, 0.0);

  x[1] = 3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -3.9999;
  EXPECT_NEAR(0., barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -4.0001;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));
  x[1] = -4.9999;
  EXPECT_NEAR(-1.234, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(barrier.inside(x));

  x[1] = 5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
  x[1] = -5.0001;
  EXPECT_NEAR(NUM_INF, barrier.energy(x), feasst::DTOL);
  EXPECT_TRUE(!barrier.inside(x));
}
