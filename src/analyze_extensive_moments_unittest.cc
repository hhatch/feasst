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
#include <limits.h>
#include "analyze_extensive_moments.h"
#include "pair_lj.h"
#include "pair_hybrid.h"
#include "mc_wltmmc.h"
#include "criteria.h"

#include "mc_wltmmc.h"
#include "ui_abbreviated.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_md.h"
#include "trial_transform.h"

using namespace feasst;

TEST(AnalyzeExtensiveMoments, constructor) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "../forcefield/data.lj"}});
  pairLJ.initData("../forcefield/data.ljb");

  AnalyzeExtensiveMoments em (&pairLJ);
  EXPECT_EQ(static_cast<int>(em.nMolTypes()), 2);

  // default values
  EXPECT_EQ (em.nFreq(), 1);
  EXPECT_EQ (em.nFreqPrint(), 1);
  EXPECT_EQ (em.safe(), false);
  EXPECT_EQ (static_cast<int>(em.order()), 0);
}

TEST(AnalyzeExtensiveMoments, setOrderParam) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";
  pairLJ.initData(molA);
  pairLJ.initData(molB);

  feasst::PairHybrid ipair(&space);
  ipair.addPair(&pairLJ);
  ipair.initEnergy();

  string orderParam = "nmol";
  double orderMin = 50, orderMax = 100, dbin = 1;
  int nbin = 51;

  const double beta = 1.0;
  feasst::CriteriaWLTMMC cc(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  cc.addActivity(exp(-1.));
  feasst::WLTMMC mc (&space, &ipair, &cc);

  AnalyzeExtensiveMoments em (&ipair);
  em.setOrderParam (3, &cc);

  EXPECT_EQ (em.safe(), true);
  EXPECT_EQ (static_cast<int>(em.order()), 3);
}

TEST(AnalyzeExtensiveMoments, writeRestart) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";
  pairLJ.initData(molA);
  pairLJ.initData(molB);

  feasst::PairHybrid ipair(&space);
  ipair.addPair(&pairLJ);
  ipair.initEnergy();

  string orderParam = "nmol";
  double orderMin = 50, orderMax = 100, dbin = 1;
  int nbin = 51;

  const double beta = 1.0;
  feasst::CriteriaWLTMMC cc(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  cc.addActivity(exp(-1.));
  feasst::WLTMMC mc (&space, &ipair, &cc);

  AnalyzeExtensiveMoments em (&ipair);
  em.setOrderParam (3, &cc);
  em.writeRestart("tmp/extmom_rst.dat");
}

TEST(AnalyzeExtensiveMoments, update1) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";

  pairLJ.initData(molA);
  pairLJ.initData(molB);

  pairLJ.linearShift(true);

  feasst::PairHybrid pair(&space);
  pair.addPair(&pairLJ);
  pair.initEnergy();

  const double beta = 1.;
  int orderMin = 50, orderMax = 60, dbin = 1;
  unsigned int nbin = 11;
  string orderParam = "nmol";
  feasst::CriteriaWLTMMC criteria(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  criteria.addActivity(exp(-1.));
  feasst::MC mc(&space, &pair, &criteria);

  AnalyzeExtensiveMoments em (&pair);
  unsigned int order = 3;

  EXPECT_FALSE (em.safe());
  em.setOrderParam (order, &criteria);
  EXPECT_TRUE (em.safe());

  EXPECT_EQ (static_cast<int>(em.nMolTypes()), 2);

  feasst::deleteTrial(&mc, molA.c_str());
  feasst::deleteTrial(&mc, molB.c_str());
    feasst::addTrial(&mc, molA.c_str());
  feasst::addTrial(&mc, molB.c_str());

  mc.nMolSeek(orderMin, molA.c_str());
  em.update();

  // everything to 0 power is 1 at first bin (ntot = 50) which has just been recorded
  auto mom = em.getMoments();
  int idx = 0;
  EXPECT_TRUE (fabs(mom[0].average() - 1) < DTOL);
  for (unsigned int i = 0; i < em.nMolTypes(); ++i) {
    for (unsigned int j = 0; j <= order; ++j) {
      for (unsigned int k = 0; k < em.nMolTypes(); ++k) {
        for (unsigned int m = 0; m <= order; ++m) {
          for (unsigned int p = 0; p <= order; ++p) {
            // other bins have not been recorded so should all return 0
            for (unsigned int n = 1; n < nbin; ++n) {
              idx = em.getIdx(i,j,k,m,p,n+orderMin);
              EXPECT_TRUE (fabs(mom[idx].average() - 0) < DTOL);
            }
          }
        }
      }
    }
  }

  mc.nMolSeek(orderMax, molB.c_str()); // fill to orderMax
  em.update();
  mom = em.getMoments();

  // everything to 0 power is 1 at first bin (ntot = 50) and last bin (ntot = 100) which has just been recorded
  idx = em.getIdx(0,0,0,0,0,orderMax);
  EXPECT_TRUE (fabs(mom[idx].average() - 1) < DTOL);
  idx = em.getIdx(0,0,0,0,0,orderMin);
  EXPECT_TRUE (fabs(mom[idx].average() - 1) < DTOL);
  for (unsigned int i = 0; i < em.nMolTypes(); ++i) {
    for (unsigned int j = 0; j <= order; ++j) {
      for (unsigned int k = 0; k < em.nMolTypes(); ++k) {
        for (unsigned int m = 0; m <= order; ++m) {
          for (unsigned int p = 0; p <= order; ++p) {
            // other bins have not been recorded so should all return 0
            for (unsigned int n = 1; n < nbin-1; ++n) {
              idx = em.getIdx(i,j,k,m,p,n+orderMin);
              EXPECT_TRUE (fabs(mom[idx].average() - 0) < DTOL);
            }
          }
        }
      }
    }
  }
}

TEST(AnalyzeExtensiveMoments, readRestart_manual1) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";

  pairLJ.initData(molA);
  pairLJ.initData(molB);

  pairLJ.linearShift(true);

  feasst::PairHybrid pair(&space);
  pair.addPair(&pairLJ);
  pair.initEnergy();

  const double beta = 1.;
  double orderMin = 50, orderMax = 100, dbin = 1;
  unsigned int nbin = 51;
  string orderParam = "nmol";
  feasst::CriteriaWLTMMC criteria(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  criteria.addActivity(exp(-1.));
  feasst::MC mc(&space, &pair, &criteria);

  AnalyzeExtensiveMoments em (&pair);
  int order = 3;
  em.setOrderParam (order, &criteria);

  feasst::deleteTrial(&mc, molA.c_str());
  feasst::deleteTrial(&mc, molB.c_str());
    feasst::addTrial(&mc, molA.c_str());
  feasst::addTrial(&mc, molB.c_str());

  mc.nMolSeek(50, molA.c_str());
  em.update();

  mc.nMolSeek(100, molB.c_str()); // add 50 MORE for 100 total
  em.update();

  // print
  string rName = "tmp/extmom2_rst.dat";
  em.writeRestart(rName.c_str());

  AnalyzeExtensiveMoments em2 (&pair);
  em2.setOrderParam (order, &criteria);

  // manually restart
  em2.readRestart(rName.c_str()); // before crit assigned
  em2.setCrit(&criteria);
}

TEST(AnalyzeExtensiveMoments, readRestart_manual2) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";

  pairLJ.initData(molA);
  pairLJ.initData(molB);

  pairLJ.linearShift(true);

  feasst::PairHybrid pair(&space);
  pair.addPair(&pairLJ);
  pair.initEnergy();

  const double beta = 1.;
  double orderMin = 50, orderMax = 100, dbin = 1;
  unsigned int nbin = 51;
  string orderParam = "nmol";
  feasst::CriteriaWLTMMC criteria(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  criteria.addActivity(exp(-1.));
  feasst::MC mc(&space, &pair, &criteria);

  AnalyzeExtensiveMoments em (&pair);
  int order = 3;
  em.setOrderParam (order, &criteria);

  feasst::deleteTrial(&mc, molA.c_str());
  feasst::deleteTrial(&mc, molB.c_str());
    feasst::addTrial(&mc, molA.c_str());
  feasst::addTrial(&mc, molB.c_str());

  mc.nMolSeek(50, molA.c_str());
  em.update();

  mc.nMolSeek(100, molB.c_str()); // add 50 MORE for 100 total
  em.update();

  // print
  string rName = "tmp/extmom2_rst.dat";
  em.writeRestart(rName.c_str());

  AnalyzeExtensiveMoments em2 (&pair);
  em2.setOrderParam (order, &criteria);

  // manually restart
  em2.setCrit(&criteria);
  em2.readRestart(rName.c_str()); // after crit assigned
}

TEST(AnalyzeExtensiveMoments, readRestart_constructor) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";

  pairLJ.initData(molA);
  pairLJ.initData(molB);

  pairLJ.linearShift(true);

  feasst::PairHybrid pair(&space);
  pair.addPair(&pairLJ);
  pair.initEnergy();

  const double beta = 1.;
  double orderMin = 50, orderMax = 100, dbin = 1;
  unsigned int nbin = 51;
  string orderParam = "nmol";
  feasst::CriteriaWLTMMC criteria(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  criteria.addActivity(exp(-1.));
  feasst::MC mc(&space, &pair, &criteria);

  AnalyzeExtensiveMoments em (&pair);
  int order = 3;
  em.setOrderParam (order, &criteria);

  feasst::deleteTrial(&mc, molA.c_str());
  feasst::deleteTrial(&mc, molB.c_str());
    feasst::addTrial(&mc, molA.c_str());
  feasst::addTrial(&mc, molB.c_str());

  mc.nMolSeek(50, molA.c_str());
  em.update();

  mc.nMolSeek(100, molB.c_str()); // add 50 MORE for 100 total
  em.update();

  // print
  string rName = "tmp/extmom2_rst.dat";
  em.writeRestart(rName.c_str());

  AnalyzeExtensiveMoments em2 (&pair, rName.c_str());
  em2.setCrit(&criteria);
}

TEST(AnalyzeExtensiveMoments, readRestart_check) {
  feasst::Space space(3);
  space.initBoxLength(20);
  feasst::PairLJ pairLJ(&space, {{"rCut", "3."}, {"molType", "none"}});

  string molA = "../forcefield/data.lj";
  string molB = "../forcefield/data.ljb";

  pairLJ.initData(molA);
  pairLJ.initData(molB);

  pairLJ.linearShift(true);

  feasst::PairHybrid pair(&space);
  pair.addPair(&pairLJ);
  pair.initEnergy();

  const double beta = 1.;
  double orderMin = 50, orderMax = 100, dbin = 1;
  unsigned int nbin = 51;
  string orderParam = "nmol";
  feasst::CriteriaWLTMMC criteria(beta, exp(-1), orderParam.c_str(), orderMin-0.5*dbin, orderMax+0.5*dbin, nbin);
  criteria.addActivity(exp(-1.));
  feasst::MC mc(&space, &pair, &criteria);

  AnalyzeExtensiveMoments em (&pair);
  unsigned int order = 3;
  em.setOrderParam (order, &criteria);

  feasst::deleteTrial(&mc, molA.c_str());
  feasst::deleteTrial(&mc, molB.c_str());
    feasst::addTrial(&mc, molA.c_str());
  feasst::addTrial(&mc, molB.c_str());

  mc.nMolSeek(50, molA.c_str());
  em.update();

  mc.nMolSeek(100, molB.c_str()); // add 50 MORE for 100 total
  em.update();

  // print
  string rName = "tmp/extmom2_rst.dat";
  em.writeRestart(rName.c_str());

  AnalyzeExtensiveMoments em2 (&pair, rName.c_str());
  em2.setCrit(&criteria);

  auto mom = em.getMoments(), mom2 = em2.getMoments();
  for (unsigned int i = 0; i < em.nMolTypes(); ++i) {
    for (unsigned int j = 0; j <= order; ++j) {
      for (unsigned int k = 0; k < em.nMolTypes(); ++k) {
        for (unsigned int m = 0; m <= order; ++m) {
          for (unsigned int p = 0; p <= order; ++p) {
            for (unsigned int n = 0; n < nbin; ++n) {
              int idx = em.getIdx(i,j,k,m,p,n+orderMin);
              int idx2 = em2.getIdx(i,j,k,m,p,n+orderMin);
              EXPECT_EQ (idx, idx2);
              EXPECT_NEAR (mom[idx].average(), mom2[idx2].average(), DTOL);
              /*EXPECT_NEAR (mom[idx].nValues(), mom2[idx].nValues(), DTOL);
              EXPECT_NEAR (mom[idx].sumSq(), mom2[idx].sumSq(), DTOL);
              EXPECT_NEAR (mom[idx].sum(), mom2[idx].sum(), DTOL);*/
            }
          }
        }
      }
    }
  }
}
