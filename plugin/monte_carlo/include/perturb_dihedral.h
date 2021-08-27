
#ifndef FEASST_MONTE_CARLO_PERTURB_DIHEDRAL_H_
#define FEASST_MONTE_CARLO_PERTURB_DIHEDRAL_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

/**
  Similar to PerturbDistanceAngle, except that the truncated conical shell has
  the additional constraint from the dihedral angle formed from the mobile site
  and its three anchors.
  Currently implemented for harmonic bonds (exponent: 2), but could add an
  optional exponent model parameter to generalize this.
 */
class PerturbDihedral : public PerturbDistanceAngle {
 public:
  explicit PerturbDihedral(argtype args = argtype());
  explicit PerturbDihedral(argtype * args);

  /// Same as PerturbDistance.
  /// If 2D, angles are positive when clockwise.
  /// Thus, when 2D and reverse (e.g., kji instead of ijk), angle = 2pi - angle.
  void precompute(TrialSelect * select, System * system) override;

  /// Return the dihedral angle.
  double dihedral() const { return dihedral_; }

  /// Return the spring constant. If rigid, return -1 (default).
  double spring_constant() const { return spring_constant_; }

  /// Return true if dihedral is rigid (e.g., if spring_constant is -1).
  bool is_rigid() const;

  /// Return the randomly selected angle from the potential.
  /// If the spring constant is -1 (rigid), simply return the angle.
  double random_dihedral(Random * random,
    const double beta,  /// inverse temperature
    const int dimension) const;

  /// Place mobile site with given bond distance, angle and dihedral.
  void place_dihedral(const double distance, const double angle,
    const double dihedral,
    System * system,
    TrialSelect * select);

  void move(System * system,
    TrialSelect * select,
    Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDihedral(std::istream& istr);
  virtual ~PerturbDihedral() {}

 protected:
  void serialize_perturb_dihedral_(std::ostream& ostr) const;

 private:
  double dihedral_ = 0.;
  double spring_constant_ = -1;

  // temporary
  Position rjk_;
  Position rkl_;
  Position origin_;
  RotationMatrix rot_mat_;
};

inline std::shared_ptr<PerturbDihedral> MakePerturbDihedral(
    argtype args = argtype()) {
  return std::make_shared<PerturbDihedral>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DIHEDRAL_H_
