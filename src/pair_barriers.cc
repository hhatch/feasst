/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_barriers.h"

namespace feasst {

PairBarriers::PairBarriers(Space* space, SpeciesBarriers* specbarriers)
  : Pair(space),
    specbarriers_(specbarriers) {
  defaultConstruction_();
}

/*!
 * Restart from file.
 *
 * \param [in] space Space class to use
 * \param [in] fileName Name of file to re-instantiate class from
 */
PairBarriers::PairBarriers (Space* space, const char* fileName) : Pair(space, fileName) {
  defaultConstruction_();

  deWall_ = fstod("deWall", fileName);
  string fname = fstos("SpecBarriersFileName", fileName);

  dummySpecBarr_ = new SpeciesBarriers (fname.c_str());
    specbarriers_ = dummySpecBarr_;
}

/**
 * defaults in constructor
 */
void PairBarriers::defaultConstruction_() {
  className_.assign("PairBarriers");
  dummySpecBarr_ = nullptr;
}

/**
 * write restart file
 */
void PairBarriers::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  string fn = fileName;
  fn.append("_specbarr");
  file << "# deWall " << deWall_ << endl;
  file << "# SpecBarriersFileName " << fn.c_str() << endl;
  specbarriers_->writeRestart(fn.c_str());
}

double PairBarriers::multiPartEner(const vector<int> mpart, const int flag) {
  if (flag == 0) {}; //remove unused parameter warning
  const vector<double> &x = space_->x();
  peSRone_ = 0;
  for (int impart = 0; impart < int(mpart.size()); ++impart) {
    const int iPart = mpart[impart];
    const int iType = space_->type()[iPart];

    vector<double> coord(x.begin() + iPart*dimen_,
                         x.begin() + (iPart+1)*dimen_);
    peSRone_ += specbarriers_->energy(iType, coord);
  }
  return peSRone_;
}

void PairBarriers::initEnergy() {
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  // compute energy
  multiPartEner(space_->listAtoms());
  peTot_  = peSRone_;
}

void PairBarriers::update(const vector<int> mpart,    //!< particles involved in move
                    const int flag,         //!< type of move
                    const char* uptype    //!< description of update type
  ) {
  std::string uptypestr(uptype);
  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deWall_ = peSRone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peTot_ += peSRone_ - deWall_;
    }
    if (flag == 2) {
      peTot_ -= deWall_;
    }
    if (flag == 3) {
      peTot_ += deWall_;
    }
  }
}

shared_ptr<PairBarriers> makePairBarriers(Space* space,
  SpeciesBarriers* specbarriers) {
  return make_shared<PairBarriers>(space, specbarriers);
}

}  // namespace feasst
