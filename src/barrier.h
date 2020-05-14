/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef BARRIER_H_
#define BARRIER_H_

#include "./base.h"
#include "./pair.h"

namespace feasst {

// Virtual class for different features of the barrier.
class Feature : public Base {
 public:
  Feature ();
  Feature (const char* fileName) { ASSERT(0, "attempting to restart a feature without an implemented restart constructor"); }
  virtual ~Feature() {;}
  virtual bool safe (const Pair& pair) = 0;
  virtual double energy (const vector < double > &coordinate) = 0;
  virtual bool inside (const vector < double > &coordinate) = 0;
  virtual double rCut () = 0;
  virtual void writeRestart (const char* fileName) = 0;

  void writeRestartBase (const char* fileName) {
    fileBackUp(fileName);
    std::ofstream file(fileName);
    file << "# className " << className_ << endl;
    file.close();
  }

  /// Convert the Square well slit pore into a Lennard-Jones slit pore.
  /// /f$ U_{wall} = \epsilon \left[ A*(\sigma/x)^{b1} - (\sigma/x)^{b2} \right]\f$
  void initLJ(const double A, const double b1, const double b2,
              /// Put the LJ interaction on only the lower wall if == 1
              const int lowerOnly = 0);

 protected:
  // Lennard-Jones parameters
  int ljFlag_ = 0, ljLowerOnly_;
  double ljA_, ljb1_, ljb2_;
};

// Overall, composite barrier composed of different features for a single particle type.
class Barrier : public Base {
public:
  Barrier (); //!< Constructor
  Barrier(const char* fileName);
  virtual ~Barrier() {}

  bool safe (const Pair& pair);
  double energy (const vector < double > &coordinate);
  bool inside (const vector < double > &coordinate);
  double rCut ();

  void writeRestart (const char* fileName);

  void addHardSlitPore (const double upper, const double lower, const unsigned int dimen);
  void addHardSlitPore (const char* fileName);

  void addSqwSlitPore (const double upper, const double lower, const double eps, const double width, const unsigned int dimen);
  void addSqwSlitPore (const char* fileName);

  void addSqwCylinder (const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen);
  void addSqwCylinder (const char* fileName);

  // add additional factory functions here ...

  vector < shared_ptr < Feature > > features() { return features_; }
  shared_ptr < Feature > features(const int index) { return features_[index]; }

protected:
  vector < shared_ptr < Feature > > features_;

  void add_ (shared_ptr < Feature > feature);
};

/**
 * This class serves as a container for derives classes with specific barrier
 * implementations, in order to interface with PairBarriers.
 */
class SpeciesBarriers : public Base {
public:
  SpeciesBarriers (const int nTypes) { set_(nTypes); } //!< Constructor
  SpeciesBarriers (const char* fileName);
  ~SpeciesBarriers() {;}

  double energy (const unsigned int type, const vector < double > &coordinate);
  bool inside (const unsigned int type, const vector < double > &coordinate);
  bool safe (const Pair& pair);
  double rCut (const unsigned int type);

  void writeRestart (const char* fileName);

  void addHardSlitPore (const unsigned int type, const double upper, const double lower, const unsigned int dimen);
  void addSqwSlitPore (const unsigned int type, const double upper, const double lower, const double eps, const double width, const unsigned int dimen);
  void addSqwCylinder (const unsigned int type, const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen);

  // add additional feature factory functions here ...

  vector < shared_ptr < Barrier > > specbarriers() { return specbarriers_; }
  shared_ptr < Barrier > specbarriers(const int index) { return specbarriers_[index]; }

protected:
  vector < shared_ptr < Barrier > > specbarriers_;
  vector < bool > assigned_;

  void set_ (const int nTypes);
  void add_ (const unsigned int type, shared_ptr < Barrier > barrier);
};

// "Hard" slit pore with no interactions with the species
class HardSlitPore : public Feature {
public:
  HardSlitPore (const char* fileName);
  HardSlitPore (const double upper, const double lower, const unsigned int dimen) { set_(upper, lower, dimen); }
  ~HardSlitPore () {;}

  void writeRestart (const char* fileName);

  bool safe (const Pair& pair);
  double energy (const vector < double > &coordinate);
  bool inside (const vector < double > &coordinate);
  double rCut () { return 0; }

  double getUpper () const { return upper_; }
  double getLower () const { return lower_; }
  unsigned int getDimen () const { return dimen_;}

 protected:
  double upper_, lower_;
  unsigned int dimen_;

  void set_ (const double upper, const double lower, const unsigned int dimen);
};

// Slit pore with square well interactions
class SqwSlitPore : public Feature {
public:
  SqwSlitPore (const char* fileName);
  SqwSlitPore (const double upper, const double lower, const double eps, const double width, const unsigned int dimen) { set_(upper, lower, eps, width, dimen);}
  ~SqwSlitPore () {;}

  void writeRestart (const char* fileName);

  bool safe (const Pair& pair);
  double energy (const vector < double > &coordinate);
  bool inside (const vector < double > &coordinate);
  double rCut () { return width_; }

  double getUpper () const { return upper_; }
        double getLower () const { return lower_; }
  double getEps () const { return eps_; }
  double getWidth () const { return width_; }
        unsigned int getDimen () const { return dimen_;}

 protected:
  double upper_, lower_, eps_, width_;
  unsigned int dimen_;

  void set_ (const double upper, const double lower, const double eps, const double width, const unsigned int dimen);
};

// Slit pore with square well interactions
class SqwCylinder : public Feature {
public:
  SqwCylinder (const char* fileName);
  SqwCylinder (const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen) { set_(pt, radius, eps, width, dimen); }
  ~SqwCylinder () {;}

  void writeRestart (const char* fileName);

  bool safe (const Pair& pair);
  double energy (const vector < double > &coordinate);
  bool inside (const vector < double > &coordinate);
  double rCut () { return width_; }

  double getRadius () const { return radius_; }
        double getEps () const { return eps_; }
  double getWidth () const { return width_; }
        unsigned int getDimen () const { return dimen_;}
  const vector < double > getPt () const { return pt_; }

 protected:
  vector < double > pt_;
  vector < double > vec_;
  double radius_, eps_, width_;
  unsigned int dimen_;

  void set_ (const vector < double > pt, const double radius, const double eps, const double width, const unsigned int dimen);
};

// Add additional barriers here ...

}  // namespace feasst

#endif  // BARRIER_H_
