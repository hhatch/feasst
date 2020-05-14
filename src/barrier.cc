/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./barrier.h"
#include <string>

namespace feasst {

/*!
 * Instantiate a species barrier from a restart file.
 *
 * \param [in] fileName Name of restart file.
 */
SpeciesBarriers::SpeciesBarriers (const char* fileName) {
  ASSERT(fileExists(fileName), "cannot locate restart file (" << fileName << ") for SpeciesBarriers");
  string name = fstos("className", fileName);
  ASSERT (name.compare("SpeciesBarriers") == 0, "file " << fileName << " does not seem to correspond to a SpeciesBarriers");

  int nAtomTypes = fstoi("nAtomTypes", fileName);
  set_(nAtomTypes);

  for (unsigned int i = 0; i < static_cast<unsigned int>(nAtomTypes); ++i) {
    string nn = "assigned";
    nn.append(std::to_string(i));
    int assigned = fstoi(nn.c_str(), fileName);
    assigned_[i] = assigned;
    if (assigned) {
      // Do not need to invoke "incremental" factory functions here since the restarted barrier should be "complete" for each species.
      string fn = "rstFileName";
      fn.append(std::to_string(i));
      string cn = fstos(fn.c_str(), fileName);
      shared_ptr < Barrier > newBarrier = make_shared < Barrier > (cn.c_str());
      specbarriers_[i] = newBarrier;
    }
  }
}

/*!
 * Instantiate a species barrier.
 *
 * \param [in] nTypes Number of atom types.
 */
void SpeciesBarriers::set_ (const int nTypes) {
  className_.assign("SpeciesBarriers");
  ASSERT (nTypes > 0, "number of particle types < 1");
  specbarriers_.resize(nTypes);
  assigned_.resize(nTypes);
  assigned_.assign(nTypes, false);
}

/*!
 * Write a restart file.
 *
 * \param [in] fileName Name of restart file.
 */
void SpeciesBarriers::writeRestart (const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# nAtomTypes " << specbarriers_.size() << endl;
  for (unsigned int i = 0; i < specbarriers_.size(); ++i) {
    file << "# assigned" << i << " " << int(assigned_[i]) << endl;
    string fn = fileName;
    fn.append("_barrier");
    fn.append(std::to_string(i));
    file << "# rstFileName" << i << " " << fn << endl;
    if (assigned_[i] == 1) {
      specbarriers_[i]->writeRestart(fn.c_str());
    }
  }
}

/*!
 * Check that the barriers present to do not violate any obvious problems with the boundaries.
 *
 * \param [in] pair Pair class containing all pairwise interactions for the system.
 *
 * \return True if safe, False if there is an error.
 */
bool SpeciesBarriers::safe (const Pair& pair) {
  bool safe = true;

  for (vector < shared_ptr < Barrier > >::iterator it = specbarriers_.begin(); it != specbarriers_.end(); ++it) {
    if (!(*it)->safe(pair)) {
      safe = false;
      break;
    }
  }

  return safe;
}

/*!
 * Get the rCut for a given atom or particle type.  Defaults ot 0 if no barriers for a given atom type.
 *
 * \param [in] type Atom type.
 *
 * \return Max r_{cut} of barrier's features with atom type.
 */
double SpeciesBarriers::rCut (const unsigned int type) {
  ASSERT (type < specbarriers_.size(), "illegal barrier type");

  if (assigned_[type]) {
    return specbarriers_[type]->rCut();
  }

  return 0.0;
}

/*!
 * Assign a barrier to the class. This is used for creating the first barrier for a given atom type.
 *
 * \param [in] type Atom type index this barrier applies to.
 * \param [in] barrier Barrier itself.
 */
void SpeciesBarriers::add_ (const unsigned int type, shared_ptr < Barrier > barrier) {
  ASSERT (type < specbarriers_.size(), "illegal barrier type");
  specbarriers_[type] = barrier;
  assigned_[type] = true;
}

/*!
 * Get the total potential energy resulting from interaction with all features in the barrier for a given atom type.
 *
 * \param [in] int Atom type.
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return pe Total potential energy.
 */
double SpeciesBarriers::energy (const unsigned int type, const vector < double > &coordinate) {
  double pe = 0.0;

  if (assigned_[type]) {
    pe = specbarriers_[type]->energy(coordinate);
  }

  return pe;
}

/*!
 * Factory function to add a "hard" slit pore feature to a species barrier.
 *
 * \param [in] type Index of atom type this applies to (starts from 0).
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] dimen Cartesian dimension the slit pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void SpeciesBarriers::addHardSlitPore (const unsigned int type, const double upper, const double lower, const unsigned int dimen) {
  if (assigned_[type]) {
    // Barrier already exists
    specbarriers_[type]->addHardSlitPore(upper, lower, dimen);
  } else {
    // Must create a new (empty) Barrier, then add to it.
    shared_ptr < Barrier > newBarrier = make_shared < Barrier > ();
    newBarrier->addHardSlitPore(upper, lower, dimen);
    add_(type, newBarrier);
  }
}

/*!
 * Factory function to add a "square well" slit pore feature to a species barrier.
 *
 * \param [in] type Index of atom type this applies to (starts from 0).
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] eps Interation energy, U = -eps if in range.
 * \param [in] width Width of interaction range next to wall, e.g., U = -eps if upper > r > upper - width.
 * \param [in] dimen Cartesian dimension the slit pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void SpeciesBarriers::addSqwSlitPore (const unsigned int type, const double upper, const double lower, const double eps, const double width, const unsigned int dimen) {
  if (assigned_[type]) {
    // Barrier already exists
    specbarriers_[type]->addSqwSlitPore(upper, lower, eps, width, dimen);
  } else {
    // Must create a new (empty) Barrier, then add to it.
    shared_ptr < Barrier > newBarrier = make_shared < Barrier > ();
    newBarrier->addSqwSlitPore(upper, lower, eps, width, dimen);
    add_(type, newBarrier);
  }
}

/*!
 * Factory function to add a "square well" cylindrical pore feature to a species barrier.
 *
 * \param [in] type Index of atom type this applies to (starts from 0).
 * \param [in] pt Any 3D coordinate which lies on the cylinder's central axis.
 * \param [in] radius Radius of the cylinder.
 * \param [in] eps Interaction energy parameter, U = -eps.
 * \param [in] width Width of interaction "band", where U = -eps from radius-width < r-center < radius.
 * \param [in] dimen Cartesian dimension the cylinder's axis lies along. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void SpeciesBarriers::addSqwCylinder (const unsigned int type, const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen) {
  if (assigned_[type]) {
    // Barrier already exists
    specbarriers_[type]->addSqwCylinder(pt, radius, eps, width, dimen);
  } else {
    // Must create a new (empty) Barrier, then add to it.
    shared_ptr < Barrier > newBarrier = make_shared < Barrier > ();
    newBarrier->addSqwCylinder(pt, radius, eps, width, dimen);
    add_(type, newBarrier);
  }
}

/*!
 * Test if coordinate is inside atleast one feature for the barrier corresponding to a given particle type.
 * By default, returns True if there are no barriers specified for this type (all coordinates are valid).
 *
 * \param [in] int Particle type.
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return inside If inside at least one feature of the barrier then True, else False.
 */
bool SpeciesBarriers::inside (const unsigned int type, const vector < double > &coordinate) {
  bool inside = true;

  if (assigned_[type]) {
    inside = specbarriers_[type]->inside(coordinate);
  }

  return inside;
}

/*!
 * Instantiate a composite barrier.
 */
Barrier::Barrier () {
  className_.assign("Barrier");
}

/*!
 * Instantiate a barrier from a restart file.
 * This creates all features associated with the barrier, so the barrier is "complete".
 *
 * \param [in] fileName Name of restart file.
 */
Barrier::Barrier (const char* fileName) {
  ASSERT(fileExists(fileName), "cannot locate restart file (" << fileName << ") for Barrier");
  string name = fstos("className", fileName);
  ASSERT (name.compare("Barrier") == 0, "file " << fileName << " does not seem to correspond to a Barrier");
  className_.assign("Barrier");

  int nFeatures = fstoi("nFeatures", fileName);
  ASSERT (nFeatures >= 0, "illegal number of features (" << std::to_string(nFeatures) << ")");

  for (unsigned int i = 0; i < static_cast<unsigned int>(nFeatures); ++i) {
    string nn = "rstFileName";
    nn.append(std::to_string(i));
    string fn = fstos(nn.c_str(), fileName);

    string cn = "featureName";
    cn.append(std::to_string(i));
    string cname = fstos(cn.c_str(), fileName);

    // re-create that feature
    if (cname == "HardSlitPore") {
      addHardSlitPore (fn.c_str());
    } else if (cname == "SqwSlitPore") {
      addSqwSlitPore (fn.c_str());
    } else if (cname == "SqwCylinder") {
      addSqwCylinder (fn.c_str());
    } else {
      ASSERT (0, "unrecognized feature type (" << cname << ")");
    }
  }
}
/*!
 * Write a restart file. Each feature of the barrier is written to a separate file.
 *
 * \param [in] fileName Name of restart file.
 */
void Barrier::writeRestart (const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# nFeatures " << features_.size() << endl;
  for (unsigned int i = 0; i < features_.size(); ++i) {
    string fn = fileName;
    fn.append("_feature");
    fn.append(std::to_string(i));
    file << "# featureName" << i << " " << features_[i]->className() << endl;
    file << "# rstFileName" << i << " " << fn << endl;
    features_[i]->writeRestart(fn.c_str());
  }
}

/*!
 * Add a feature to the barrier.
 *
 * \param [in] feature Shared pointer to the feature to add.
 */
void Barrier::add_ (shared_ptr < Feature > feature) {
  features_.push_back(feature);
}

/*!
 * Get the total potential energy resulting from interaction with all features in the barrier.
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return pe Total potential energy.
 */
double Barrier::energy (const vector < double > &coordinate) {
  double pe = 0.0;

  for (vector < shared_ptr < Feature > >::iterator it = features_.begin(); it != features_.end(); ++it) {
        pe += (*it)->energy(coordinate);
    }

  return pe;
}

/*!
 * Test if coordinate is inside atleast one feature.
 *
 * \param [in] int Atom type.
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return inside If inside at least one feature of the barrier then True, else False.
 */
bool Barrier::inside (const vector < double > &coordinate) {
  for (vector < shared_ptr < Feature > >::iterator it = features_.begin(); it != features_.end(); ++it) {
    if ((*it)->inside(coordinate)) {
      return true;
    }
    }

  return false;
}

/*!
 * Get max rCut for all features present.
 *
 * \return Max r_{cut} for all features.
 */
double Barrier::rCut () {
  double rc = 0.0, rcn;
  for (vector < shared_ptr < Feature > >::iterator it = features_.begin(); it != features_.end(); ++it) {
    rcn = (*it)->rCut();
    if (rcn > rc) {
      rc = rcn;
    }
  }

  return rc;
}

/*!
 * Check that the features present to do not violate any obvious problems with the boundaries.
 *
 * \param [in] pair Pair class containing all pairwise interactions for the system.
 * \return True if safe, False if there is an error.
 */
bool Barrier::safe (const Pair& pair) {
  for (vector < shared_ptr < Feature > >::iterator it = features_.begin(); it != features_.end(); ++it) {
    if (!(*it)->safe(pair)) {
      return false;
    }
  }

  return true;
}

/*!
 * Factory function to add a "hard" slit pore feature to a composite barrier.
 *
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] dimen Cartesian dimension the slit pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void Barrier::addHardSlitPore (const double upper, const double lower, const unsigned int dimen) {
  shared_ptr < Feature > newBarrier = make_shared < HardSlitPore > (upper, lower, dimen);
  add_ (newBarrier);
}

/*!
 * Factory function to add a "hard" slit pore feature to a composite barrier from restart file.
 *
 * \param [in] fileName File containing restart information.
 */
void Barrier::addHardSlitPore (const char* fileName) {
  shared_ptr < Feature > newBarrier = make_shared < HardSlitPore > (fileName);
  add_ (newBarrier);
}

/*!
 * Factory function to add a "square well" slit pore feature to a composite barrier.
 *
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] eps Interation energy, U = -eps if in range.
 * \param [in] width Width of interaction range next to wall, e.g., U = -eps if upper > r > upper - width.
 * \param [in] dimen Cartesian dimension the slit pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void Barrier::addSqwSlitPore (const double upper, const double lower, const double eps, const double width, const unsigned int dimen) {
  shared_ptr < Feature > newBarrier = make_shared < SqwSlitPore > (upper, lower, eps, width, dimen);
  add_ (newBarrier);
}

/*!
 * Factory function to add a "square well" slit pore feature to a composite barrier from restart file.
 *
 * \param [in] fileName File containing restart information.
 */
void Barrier::addSqwSlitPore (const char* fileName) {
  shared_ptr < Feature > newBarrier = make_shared < SqwSlitPore > (fileName);
  add_ (newBarrier);
}

/*!
 * Factory function to add a "square well" cylindrical pore feature to a composite barrier.
 *
 * \param [in] pt Any 3D coordinate which lies on the cylinder's central axis.
 * \param [in] radius Radius of the cylinder.
 * \param [in] eps Interaction energy parameter, U = -eps.
 * \param [in] width Width of interaction "band", where U = -eps from radius-width < r-center < radius.
 * \param [in] dimen Cartesian dimension the cylinder's axis lies along. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void Barrier::addSqwCylinder (const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen) {
  shared_ptr < Feature > newBarrier = make_shared < SqwCylinder > (pt, radius, eps, width, dimen);
  add_ (newBarrier);
}

/*!
 * Factory function to add a "square well" cylindrical pore feature to a composite barrier from restart file.
 *
 * \param [in] fileName File containing restart information.
 */
void Barrier::addSqwCylinder (const char* fileName) {
  shared_ptr < Feature > newBarrier = make_shared < SqwCylinder > (fileName);
  add_ (newBarrier);
}

/*!
 * Instantiate a feature.
 */
Feature::Feature () {
  className_.assign("Feature");
}

/*!
* Instantiate a species barrier from a restart file.
*
* \param [in] fileName Name of restart file.
 */
HardSlitPore::HardSlitPore (const char* fileName) {
  ASSERT(fileExists(fileName), "cannot locate restart file (" << fileName << ") for HardSlitPore");
  string name = fstos("className", fileName);
  ASSERT (name.compare("HardSlitPore") == 0, "file " << fileName << " does not seem to correspond to a HardSlitPore");

  double upper = 0, lower = 0;
  unsigned int dimen = -1;

  upper = fstod("upper", fileName);
  lower = fstod("lower", fileName);
  dimen = fstoi("dimen", fileName);

  set_ (upper, lower, dimen);
}

/*!
 * Assign variables to a "hard" slit pore feature.
 *
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] dimen Cartesian dimension the slit pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void HardSlitPore::set_ (const double upper, const double lower, const unsigned int dimen) {
  ASSERT (upper > lower, "upper bound must be greater than lower bound for SlitPore");
  ASSERT (dimen < 3 && dimen >= 0, "dimension can only be: 0 (x), 1 (y), or 2 (z)")
  className_.assign("HardSlitPore");
  upper_ = upper;
  lower_ = lower;
  dimen_ = dimen;
}

/*!
* Write a restart file.
*
* \param [in] fileName Name of restart file.
 */
void HardSlitPore::writeRestart (const char* fileName) {
  writeRestartBase (fileName);
  ofstream file;
  file.open(fileName, std::ios_base::app); // append information after writeRestartBase

  file << "# upper " << upper_ << endl;
  file << "# lower " << lower_ << endl;
  file << "# dimen " << dimen_ << endl;

  file.close();
}

/*!
 * Get the total potential energy resulting from interaction with the "hard" slit pore.
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return pe Total potential energy (0 or infinity).
 */
double HardSlitPore::energy (const vector < double > &coordinate) {
  ASSERT(coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");

  if (inside(coordinate)) {
    return 0.0;
  }

  return NUM_INF;
}

/*!
 * Test if a position is "inside" the "hard" slit pore (between the confining walls).
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return inside True if inside, else False.
 */
bool HardSlitPore::inside (const vector < double > &coordinate) {
  ASSERT(coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");

  if (coordinate[dimen_] > lower_ && coordinate[dimen_] < upper_) {
    return true;
  }

  return false;
}

/*!
 * Check that the barrier is wholly inside the box and prevents atoms from interacting with their periodic images across any non-periodic dimension.
 *
 * \param [in] pair Pair class containing all pairwise interactions for the system.
 *
 * \return True if safe, False if there is an error.
 */
bool HardSlitPore::safe (const Pair &pair) {
  // Check boundaries inside the box
  if (!(upper_ < pair.space()->boxLength()[dimen_]/2.0)) return false;
  if (!(lower_ > -pair.space()->boxLength()[dimen_]/2.0)) return false;

  // Check that species do not interact across PBC
  if (!(pair.space()->boxLength()[dimen_] - (upper_ - lower_) > pair.rCutMaxAll())) return false;

  return true;
}


/*!
* Instantiate a species barrier from a restart file.
*
* \param [in] fileName Name of restart file.
 */
SqwSlitPore::SqwSlitPore (const char* fileName) {
  ASSERT(fileExists(fileName), "cannot locate restart file (" << fileName << ") for SqwSlitPore");
  string name = fstos("className", fileName);
  ASSERT (name.compare("SqwSlitPore") == 0, "file " << fileName << " does not seem to correspond to a SqwSlitPore");

  double upper = 0, lower = 0, eps = 0, width = 0;
  unsigned int dimen = -1;

  upper = fstod("upper", fileName);
  lower = fstod("lower", fileName);
  eps = fstod("eps", fileName);
  width = fstod("width", fileName);
  dimen = fstoi("dimen", fileName);

  set_ (upper, lower, eps, width, dimen);
}

/*!
 * Assign variables to a "square well" slit pore feature.
 *
 * \param [in] upper Maximum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] lower Minimum value an interacting "atom's" CENTER OF MASS is allowed.
 * \param [in] eps Interaction energy parameter, U = -eps.
 * \param [in] width Width of interaction "band", where U = -eps from upper (lower) to upper -(+) width.
 * \param [in] dimen Cartesian dimension the square well pore is in. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void SqwSlitPore::set_ (const double upper, const double lower, const double eps, const double width, const unsigned int dimen) {
  ASSERT (width >= 0, "width must be >= 0");
  ASSERT (upper > lower, "upper bound must be greater than lower bound for SqwSlitPore");
  ASSERT (dimen < 3 && dimen >= 0, "dimension can only be: 0 (x), 1 (y), or 2 (z)");
  ASSERT (upper-lower >= 2*width, "slit pore interaction range overlaps between opposing boundaries, this is not allowed");
  ASSERT (eps >= 0, "interaction energy must be >= 0");
  className_.assign("SqwSlitPore");
  upper_ = upper;
  lower_ = lower;
  eps_ = eps;
  width_ = width;
  dimen_ = dimen;
}

/*!
* Write a restart file.
*
* \param [in] fileName Name of restart file.
 */
void SqwSlitPore::writeRestart (const char* fileName) {
  writeRestartBase (fileName);
  ofstream file;
  file.open(fileName, std::ios_base::app); // append information after writeRestartBase

  file << "# upper " << upper_ << endl;
  file << "# lower " << lower_ << endl;
  file << "# eps " << eps_ << endl;
  file << "# width " << width_ << endl;
  file << "# dimen " << dimen_ << endl;

  file.close();
}

void Feature::initLJ(const double A, const double b1, const double b2,
  const int lowerOnly) {
  ljFlag_ = 1;
  ljA_ =  A;
  ljb1_ = b1;
  ljb2_ = b2;
  ljLowerOnly_ = lowerOnly;
}

/*!
 * Get the total potential energy resulting from interaction with the "square well" slit pore.
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return pe Total interaction potential energy.
 */
double SqwSlitPore::energy (const vector < double > &coordinate) {
  ASSERT(coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");

  double U = 0.0;

  if (!inside(coordinate)) {
    // out of bounds
    U = NUM_INF;
  } else if (ljFlag_ == 1) {
    U = 0.;

    // upper wall
    if (ljLowerOnly_ == 0) {
      const double pos = upper_ - coordinate[dimen_];
      const double sigr = width_/pos;
      U = eps_*(ljA_*pow(sigr, ljb1_) - pow(sigr, ljb2_));
    }

    // lower wall
    const double pos = coordinate[dimen_] - lower_;
    const double sigr = width_/pos;
    U += eps_*(ljA_*pow(sigr, ljb1_) - pow(sigr, ljb2_));

  } else if (coordinate[dimen_] > upper_ - width_) {
    // in bounds, in "top" range
    U = -eps_;
  } else if (coordinate[dimen_] < lower_ + width_) {
    // in bounds, in "bottom" range
    U = -eps_;
  } else {
    U = 0.0;
  }

  return U;
}

/*!
 * Test if a position is "inside" the "square well" slit pore (between the confining walls).
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return inside True if inside, else False.
 */
bool SqwSlitPore::inside (const vector < double > &coordinate) {
  ASSERT(coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");

  if (coordinate[dimen_] > lower_ && coordinate[dimen_] < upper_) {
    return true;
  }

  return false;
}

/*!
 * Check that the barrier is wholly inside the box and prevents atoms from interacting with their periodic images across any non-periodic dimension.
 *
 * \param [in] pair Pair class containing all pairwise interactions for the system.
 *
 * \return True if safe, False if there is an error.
 */
bool SqwSlitPore::safe (const Pair &pair) {
  // Check boundaries inside the box
  if (!(upper_ < pair.space()->boxLength()[dimen_]/2.0)) return false;
  if (!(lower_ > -pair.space()->boxLength()[dimen_]/2.0)) return false;

  // Check that species do not interact across PBC
  if (!(pair.space()->boxLength()[dimen_] - (upper_ - lower_) > pair.rCutMaxAll())) return false;

  return true;
}

/*!
* Instantiate a species barrier from a restart file.
*
* \param [in] fileName Name of restart file.
 */
SqwCylinder::SqwCylinder (const char* fileName) {
  ASSERT(fileExists(fileName), "cannot locate restart file (" << fileName << ") for SqwCylinder");
  string name = fstos("className", fileName);
  ASSERT (name.compare("SqwCylinder") == 0, "file " << fileName << " does not seem to correspond to a SqwCylinder");

  vector < double > pt (3,0);
  double radius = 0, eps = 0, width = 0;
  unsigned int dimen = -1;

  pt[0] = fstod("pt0", fileName);
  pt[1] = fstod("pt1", fileName);
  pt[2] = fstod("pt2", fileName);
  radius = fstod("radius", fileName);
  eps = fstod("eps", fileName);
  width = fstod("width", fileName);
  dimen = fstoi("dimen", fileName);

  set_ (pt, radius, eps, width, dimen);
}

/*!
 * Assign variables to a "square well" cylindrical pore feature. Cylinder's axis is always oriented along one of the dimensions.  Only use for 3D systems because this requires cross products to evaluate in the current implementation.
 *
 * \param [in] pt Any 3D coordinate which lies on the cylinder's central axis.
 * \param [in] radius Radius of the cylinder.
 * \param [in] eps Interaction energy parameter, U = -eps.
 * \param [in] width Width of interaction "band", where U = -eps from radius-width < r-center < radius.
 * \param [in] dimen Cartesian dimension the cylinder's axis lies along. Allowable values are: 0 (x), 1 (y), 2 (z).
 */
void SqwCylinder::set_ (const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen) {
  ASSERT (pt.size() == 3, "must provide 3D coordinate for cylinder axis")
  ASSERT (radius > 0, "radius must be > 0");
  ASSERT (width >= 0, "width must be >= 0");
  ASSERT (radius >= width, "radius must be larger than or equal to the interaction width");
  ASSERT (dimen < 3 && dimen >= 0, "dimension can only be: 0 (x), 1 (y), or 2 (z)");
  ASSERT (eps >= 0, "interaction energy must be >= 0");

  pt_.resize(pt.size(), 0);
  pt_ = pt;
  radius_ = radius;
  eps_ = eps;
  width_ = width;
  dimen_ = dimen;

  vec_.resize(pt.size(), 0);
  vec_[dimen] = 1;
  className_.assign("SqwCylinder");
}

/*!
* Write a restart file.
*
* \param [in] fileName Name of restart file.
 */
void SqwCylinder::writeRestart (const char* fileName) {
  writeRestartBase (fileName);
  ofstream file;
  file.open(fileName, std::ios_base::app); // append information after writeRestartBase

  file << "# pt0 " << pt_[0] << endl;
  file << "# pt1 " << pt_[1] << endl;
  file << "# pt2 " << pt_[2] << endl;
  file << "# radius " << radius_ << endl;
  file << "# eps " << eps_ << endl;
  file << "# width " << width_ << endl;
  file << "# dimen " << dimen_ << endl;

  file.close();
}

/*!
 * Get the total potential energy resulting from interaction with the "square well" cylindrical pore.
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return pe Total interaction potential energy.
 */
double SqwCylinder::energy (const vector < double > &coordinate) {
  ASSERT (coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");
  ASSERT (coordinate.size() == 3, "SqwCylinder only valid for 3D systems");

  double U = 0.0;

  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  vector < double > dx (pt_.size(), 0);
  for (unsigned int i = 0; i < dx.size(); ++i) {
    dx[i] = pt_[i] - coordinate[i];
  }
  vector < double > a = crossProd(vec_, dx);

  double b = 0, c = 0;
  for (unsigned int i = 0; i < dx.size(); ++i) {
    b += a[i]*a[i];
    c += vec_[i]*vec_[i];
  }

  double r = sqrt(b/c);
  if (r < radius_ - width_) {
    U = 0.0;
  } else if (r < radius_) {
    U = -eps_;
  } else {
    return NUM_INF;
  }

  return U;
}

/*!
 * Test if a position is "inside" the "square well" cylindrical pore (between the confining walls).
 *
 * \param [in] coordinate Vector of cartesian coordinate to test.
 *
 * \return inside True if inside, else False.
 */
bool SqwCylinder::inside (const vector < double > &coordinate) {
  ASSERT (coordinate.size() > dimen_, "the requested dimension is not present in the coordinates provided to the potential function");
  ASSERT (coordinate.size() == 3, "SqwCylinder only valid for 3D systems");

  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  vector < double > dx (pt_.size(), 0);
  for (unsigned int i = 0; i < dx.size(); ++i) {
    dx[i] = pt_[i] - coordinate[i];
  }
  vector < double > a = crossProd(vec_, dx);

  double b = 0, c = 0;
  for (unsigned int i = 0; i < dx.size(); ++i) {
    b += a[i]*a[i];
    c += vec_[i]*vec_[i];
  }

  if (b/c < radius_*radius_) {
    return true;
  }

  return false;
}

/*!
 * Check that the barrier is wholly inside the box and prevents atoms from interacting with their periodic images across any non-periodic dimension.
 *
 * \param [in] pair Pair class containing all pairwise interactions for the system.
 *
 * \return True if safe, False if there is an error.
 */
bool SqwCylinder::safe (const Pair &pair) {
  for (unsigned int i = 0; i < pt_.size(); ++i) {
    if (i != dimen_) {
      // Check the pore's cross-section within the box
      if (!(pt_[i]+radius_ < pair.space()->boxLength()[i]/2.0) ||
          !(pt_[i]-radius_ > -pair.space()->boxLength()[i]/2.0)) {
        return false;
      }

      // Check that species do not interact across PBC
      if (!(pair.space()->boxLength()[i] - (2*radius_) > pair.rCutMaxAll())) return false;
    }
  }

  return true;
}

}  // namespace feasst
