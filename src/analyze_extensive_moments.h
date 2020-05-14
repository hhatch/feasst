/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

/**
 * \file analyze_extensive_moments.h
 *
 * \brief Store the extensive moments along the sampling order parameter of the system.  These are used to perform histogram extrapolation.
 *
 */

#ifndef ANALYZEEXTMOMENTS_H_
#define ANALYZEEXTMOMENTS_H_

#include "./analyze.h"
#include "./criteria_wltmmc.h"
#include <iomanip>

namespace feasst {

#define EXTPRECISION 18 // Number of digits of precision to print extensive moments out as

class AnalyzeExtensiveMoments : public Analyze {
public:
    AnalyzeExtensiveMoments (Pair *pair) : Analyze(pair) {
        defaults_ (space());
    }

    AnalyzeExtensiveMoments ( Pair *pair, const char* fileName) : Analyze (pair, fileName) {
        defaults_ (space());
        readRestart (fileName);
    }

  shared_ptr <AnalyzeExtensiveMoments> cloneShrPtr(Pair* pair) const {
    return (std::static_pointer_cast<AnalyzeExtensiveMoments, Analyze>(cloneImpl(pair)));
  }

  AnalyzeExtensiveMoments* clone (Pair* pair) const {
    AnalyzeExtensiveMoments* a = new AnalyzeExtensiveMoments(*this);
    a->reconstruct(pair); return a;
  }

  ~AnalyzeExtensiveMoments() { if (cRestart_ != nullptr) delete cRestart_; }

    void setOrderParam (const unsigned int order, CriteriaWLTMMC* crit);
  void setCrit (CriteriaWLTMMC* crit);
  void setOrder (const unsigned int order);

  void update ();
  void update (const int iMacro);

  void writeRestart (const char* fileName);
  void readRestart (const char* fileName);

  void write ();
  void write (CriteriaWLTMMC *c);

  unsigned int nMolTypes () const { return nspec_; }
  unsigned int order () const { return order_; }
  unsigned int getIdx (const int i, const int j, const int k, const int m, const int p, const int N);
  bool safe () const { return critAssigned_; }

    vector < Accumulator > getMoments () const { return moments_; }
  Accumulator getMoment (const unsigned int i) const { return moments_[i]; }

private:
    unsigned int order_, nspec_;
    bool critAssigned_;
    CriteriaWLTMMC* crit_; //!< Pointer to criteria used to sample
  CriteriaWLTMMC* cRestart_; //!< Criteria used if restarted
    vector < int > widths_;
    vector < Accumulator > moments_;

  void defaults_ (Space *space);
    void init_ (const unsigned int nBin);

  virtual shared_ptr<Analyze> cloneImpl(Pair *pair) const {
    shared_ptr<AnalyzeExtensiveMoments> a = make_shared<AnalyzeExtensiveMoments>(*this);
    a->reconstruct(pair); return a;
  }
};

/// Factory method
shared_ptr<AnalyzeExtensiveMoments> makeAnalyzeExtensiveMoments(Pair *pair);

}  // namespace feasst

#endif // ANALYZEEXTMOMENTS_H_
