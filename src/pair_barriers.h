/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_BARRIERS_H_
#define PAIR_BARRIERS_H_

#include "./pair.h"
#include "./barrier.h"

namespace feasst {

/**
 * Compute interactions of particles with walls defined by the Barrier class.
 */
class PairBarriers : public Pair {
 public:
  /// Constructor.
  PairBarriers (Space* space, SpeciesBarriers* specbarriers);

  /*
   * \return potential energy of a selection of particles multiPart, and store
   * this the Pair class variable peSRone
   */
  double multiPartEner (const vector < int > multiPart, const int flag = 0);

  /*
   * Compute the potential energy of all particles and store this in the Pair
   * class variable peTot
   */
  void initEnergy ();

  // stores, restores or updates variables to avoid order recompute
  //   of entire configuration after every change
  virtual void update (const vector < int > mpart, const int flag,
                      const char* uptype);

  PairBarriers (Space* space, const char* fileName);
  ~PairBarriers () { if (dummySpecBarr_ != nullptr) delete dummySpecBarr_; }
  virtual PairBarriers* clone (Space* space) const {
    PairBarriers* p = new PairBarriers (*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  SpeciesBarriers* specbarriers_;
  SpeciesBarriers* dummySpecBarr_;
  double deWall_ = 0.;

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairBarriers> makePairBarriers (Space* space,
  SpeciesBarriers* specbarriers);

}  // namespace feasst

#endif  // PAIR_BARRIERS_H_
