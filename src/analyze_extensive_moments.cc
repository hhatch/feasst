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
 * \file analyze_extensive_moments.cc
 *
 * \brief Measure the moments of extensive quantities during a WLTMMC simulation for extrapolation later.
 */

#include "./analyze_extensive_moments.h"

namespace feasst {

/*!
 * Set default values.
 *
 * \param [in] space Space used during simulation.
 */
void AnalyzeExtensiveMoments::defaults_ (Space *space) {
    critAssigned_ = false;
    crit_ = nullptr;
  cRestart_ = nullptr;
    nspec_ = space->nMolTypes();
    order_ = 0;
  className_ = "AnalyzeExtensiveMoments";
    ASSERT (nspec_ > 0, "number of molecule types must be > 0");
}

/*!
 * Assign the sampling order parameter information.
 *
 * \param [in] order Maximum exponent (moment) to record up to.
 * \param [in] crit CriteriaWLTMMC used for sampling.
 */
void AnalyzeExtensiveMoments::setOrderParam (const unsigned int order, CriteriaWLTMMC* crit) {
  setOrder (order);
    setCrit (crit);
}

/*!
 * Set sampling criteria.
 *
 * \param [in] crit CriteriaWLTMMC used for sampling.
 */
void AnalyzeExtensiveMoments::setCrit (CriteriaWLTMMC* crit) {
  ASSERT (crit->mType().compare("nmol") == 0, "extensive moments analysis only valid when the sampling order parameter is nmol"); // later will extend to N1 also

  // in case already restarted
  if (moments_.begin() != moments_.end()) {
        ASSERT (moments_.size() == nspec_*nspec_*(order_+1)*(order_+1)*(order_+1)*crit->nBin(), "macrostate variable's number of bins for sampling is inconsistent with restart file");
    }

  crit_ = crit;
    critAssigned_ = true;
}

/*!
 * Set maximum order (exponent) of moments to measure.
 *
 * \param [in] order Maximum exponent (moment) to record up to.
 */
void AnalyzeExtensiveMoments::setOrder (const unsigned int order) {
  ASSERT (order > 0, "order must be > 0");
  ASSERT (!(order_ > 0 and order != order_), "attempting to reassign order to a different value after it has already been set");
  order_ = order;
}

/*!
 * Record information.
 */
void AnalyzeExtensiveMoments::update() {
    ASSERT (critAssigned_, "have not assigned order parameter, see setOrderParam()");
  update(crit_->bin(space()->nMol())); // Use Ntot as the order parameter
}

/*!
 *  Get the "unrolled" index corresponding to a moment of the form N_i^j*N_k^m*U^p(N)
 *
 * \param [in] i Index of first component
 * \param [in] j Exponent on first N
 * \param [in] k Index of second component
 * \param [in] m Exponent on second N
 * \param [in] p Exponent on U
 * \param [in] N Number of molecules
 * \return idx
 */
unsigned int AnalyzeExtensiveMoments::getIdx (const int i, const int j, const int k, const int m, const int p, const int N) {
  ASSERT (critAssigned_, "have not assigned order parameter, see setOrderParam()");
  if (widths_.begin() == widths_.end()) init_ (crit_->nBin());
  return i*widths_[0] + j*widths_[1] + k*widths_[2] + m*widths_[3] + p*widths_[4] + crit_->bin(N)*widths_[5];
}

/*!
 * Update the moments with the current state of the system.
 *
 * \param [in] iMacro Value of the macrostate (Ntot)
 */
void AnalyzeExtensiveMoments::update (const int iMacro) {
  ASSERT (critAssigned_, "have not assigned order parameter, see setOrderParam()");

  // initialize the first time called
    if (moments_.begin() == moments_.end()) {
        init_ (crit_->nBin());
    }

    const vector < int > nMol = space()->nMolType();
    const long double u = pair_->peTot();
    long long unsigned int coords = 0;

    // pre-compute powers so only needs to be done once
  vector < double > n (order_+1, 1.0), u_p (order_+1, 1.0);
  vector < vector < double > > n_i (nspec_, n);
  for (unsigned int j = 1; j <= order_; ++j) {
    for (unsigned int i = 0; i < nspec_; ++i) {
      n_i[i][j] = nMol[i]*n_i[i][j-1];
    }
    u_p[j] = u*u_p[j-1];
  }

    // assign the moments to their "unrolled" position
    for (unsigned int p = 0; p <= order_; ++p) {
        for (unsigned int m = 0; m <= order_; ++m) {
            for (unsigned int k = 0; k < nspec_; ++k) {
                for (unsigned int j = 0; j <= order_; ++j) {
                    for (unsigned int i = 0; i < nspec_; ++i) {
                        coords = i*widths_[0] + j*widths_[1] + k*widths_[2] + m*widths_[3] + p*widths_[4] + (iMacro)*widths_[5];
            moments_[coords].accumulate(n_i[i][j]*n_i[k][m]*u_p[p]);
                    }
                }
            }
        }
    }
}

/*!
 * Print the moments at given increments. Prints a file to "extMom" if filename has not been set internally with initFileName().
 * This simply calls writeRestart() so printed files are the same format as restarts.
 */
void AnalyzeExtensiveMoments::write() {
  ASSERT (critAssigned_, "CriteriaWLTMMC has not been assigned to AnalyzeExtensiveMoments, cannot write")
  write (crit_);
}

/*!
 * Print the moments at given increments. Prints a file to "extMom" if filename has not been set internally with initFileName().
 * This simply calls writeRestart() so printed files are the same format as restarts.
 *
 * \param [in] c WLTMMC criteria being used to sample with - however this is disregarded and the one assigned to this class internally is used instead.
 */
void AnalyzeExtensiveMoments::write(CriteriaWLTMMC *c) {
  if (!fileName_.empty()) {
    writeRestart (fileName_.c_str());
  } else {
    writeRestart ("extMom");
  }
}

/*!
 * Write restart file to pick up from later.
 *
 * \param [in] fileName Name of restart file.
 */
void AnalyzeExtensiveMoments::writeRestart (const char* fileName) {
  writeRestartBase (fileName);

    ofstream file;
    file.open(fileName, std::ios_base::app); // append information after writeRestartBase

    file << "# orderParam " << crit_->mType() << endl;
    file << "# mMin " << std::setprecision(EXTPRECISION) << crit_->mMin() << endl;
    file << "# mMax " << std::setprecision(EXTPRECISION) << crit_->mMax() << endl;
    file << "# mBin " << std::setprecision(EXTPRECISION) << crit_->mBin() << endl;
    file << "# nBin " << crit_->nBin() << endl;
    file << "# maxOrder " << order_ << endl;
    file << "# nSpec " << nspec_ << endl;
  vector < double > l = space()->boxLength();
  double V = 1.0;
  for (vector < double >::iterator it = l.begin(); it != l.end(); ++it) {
    V *= (*it);
  }
  file << "# volume " << V << endl;
    file << "# opIdx\tnValues\tSum\tSumSq\t<N_i\t^j*\tN_k\t^m*\tU^p>" << endl;
    file << "# BEGIN" << endl;

    long long unsigned int ctr = 0;
    if (moments_.size() == nspec_*nspec_*(order_+1)*(order_+1)*(order_+1)*crit_->nBin()) {
        for (unsigned int nidx = 0; nidx < static_cast < unsigned int >(crit_->nBin()); ++nidx) {
            for (unsigned int p = 0; p <= order_; ++p) {
                for (unsigned int m = 0; m <= order_; ++m) {
                    for (unsigned int k = 0; k < nspec_; ++k) {
                        for (unsigned int j = 0; j <= order_; ++j) {
                            for (unsigned int i = 0; i < nspec_; ++i) {
                                file << nidx << "\t" << moments_[ctr].nValues() << "\t" << moments_[ctr].sum() << "\t" << moments_[ctr].sumSq() << "\t" << i << "\t" << j << "\t" << k << "\t" << m << "\t" << p << endl;
                                ctr++;
                            }
                        }
                    }
                }
            }
        }
    }

    file << "# END" << endl;
    file.close();

  // also write the criteria used
  if (critAssigned_) {
    string critFilename = fileName;
    critFilename.append("_crit");
    crit_->writeRestart(critFilename.c_str());
  }
}

/*!
 * Read restart file and re-initialize the class to pick up from wherever it was left off. The user must still call setCrit() afterwards before running a simulation.
 *
 * \param [in] fileName Name of restart file.
 */
void AnalyzeExtensiveMoments::readRestart (const char* fileName) {
    string op = fstos("orderParam", fileName);
    ASSERT (op.compare("nmol") == 0, "extensive moments analysis only valid when the sampling order parameter is nmol");
    const unsigned int nBin = fstoi("nBin", fileName);

    if (order_ > 0) {
        // was specifed in setOrderParam()
        const unsigned int x = fstoi("maxOrder", fileName);
        ASSERT (x == order_, "maximum order specified is inconsistent with restart file");
    } else {
        // if not specified already
        order_ = fstoi("maxOrder", fileName);
    }
    ASSERT (order_ > 0, "order must be > 0");
    nspec_ = fstoi("nSpec", fileName);
    ASSERT (nspec_ > 0, "number of molecule types must be > 0");
    ASSERT (static_cast < unsigned int >(space()->nMolTypes()) == nspec_, "number of molecule types specified in space inconsistent with restart file");

  init_ (nBin); // must be called after nspec_ and order_ assigned

    // find beginning and ending
    std::ifstream file(fileName);
    string line;
    string startingPt("BEGIN"), endingPt("END");
    long double endIdx = 0.0, startIdx = 0.0;
    getline(file, line);
  endIdx += 1.0;
    while ( (line.find(endingPt) == std::string::npos) && (file) ) {
      getline(file, line);
      endIdx += 1.0;
    }
    file.close();

    file.open(fileName, std::ios_base::in);
    getline(file, line);
  startIdx += 1.0;
    while ( (line.find(startingPt) == std::string::npos) && (file) ) {
      getline(file, line);
      startIdx += 1.0;
    }

    // check file size is correct
    const long long unsigned int nM = nspec_*nspec_*(order_+1)*(order_+1)*(order_+1);
    ASSERT (std::abs(endIdx - startIdx - 2 + 1 - nBin*nM) < DTOL, "restart file constains the incorrect number of lines");

    // check moments size is correct
    ASSERT (moments_.size() == nBin*nM, "bad initialization of moments");

    long unsigned int dummyIdx1;
  long double nV, nSum, nSumSq;

    // read from file
  vector < int > expo (5,0);
  int idx = 0;
    for (unsigned long long int i = 0; i < moments_.size(); ++i) {
        file >> dummyIdx1 >> nV >> nSum >> nSumSq;
        Accumulator newMom (nV, nSum, nSumSq);
        for (unsigned int mm = 0; mm < expo.size(); ++mm) {
            file >> expo[mm];
        }
    idx = expo[0]*widths_[0] + expo[1]*widths_[1] + expo[2]*widths_[2] + expo[3]*widths_[3] + expo[4]*widths_[4] + (dummyIdx1)*widths_[5]; // must do this calculation manually here since getIdx requires crit to be set, which is not the case if instantiating directly from constructor
        moments_[idx] = newMom;
    }
    file.close();

  // re-create the criteria used internally
  string critFilename = fileName;
  critFilename.append("_crit");
  cRestart_ = new CriteriaWLTMMC (critFilename.c_str());
  setCrit(cRestart_);
}

/*!
 * Initialize the vectors to store information.
 *
 * \param [in] nBin Number of bins the macrostate distribution is going to be given.
 */
void AnalyzeExtensiveMoments::init_ (const unsigned int nBin) {
    moments_.erase(moments_.begin(), moments_.end());
    moments_.resize(nBin*nspec_*nspec_*(order_+1)*(order_+1)*(order_+1));

    widths_.resize(6,0);

    // "width" of each unrolled dimension
    vector < int > nbn (6,1);
    nbn[0] = nspec_;
  nbn[1] = order_+1;
  nbn[2] = nspec_;
  nbn[3] = order_+1;
  nbn[4] = order_+1;
  nbn[5] = nBin;
    for (unsigned int i = 0; i < widths_.size(); ++i) {
        if (i == 0) {
            widths_[i] = 1;
        } else {
            widths_[i] = widths_[i-1]*nbn[i-1];
        }
    }
}

shared_ptr<AnalyzeExtensiveMoments> makeAnalyzeExtensiveMoments(Pair *pair) {
  return make_shared<AnalyzeExtensiveMoments>(pair);
}

}  // namespace feasst
