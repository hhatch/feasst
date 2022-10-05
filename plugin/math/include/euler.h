
#ifndef FEASST_MATH_EULER_H_
#define FEASST_MATH_EULER_H_

#include <string>
#include <vector>
#include "math/include/matrix.h"

namespace feasst {

/**
  The x-convention is defined as follows:
  1. rotate by phi [-pi, pi] about the z-axis.
  2. rotate by theta [0, pi] about the new x-axis.
  3. rotate by psi [-pi, pi] about the new z-axis.
  See https://mathworld.wolfram.com/EulerAngles.html.
 */
class Euler {
 public:
  Euler() {}

  /// Construct from values.
  Euler(const double phi, const double theta, const double psi) {
    set(phi, theta, psi); }

  /// Set the Euler angles.
  void set(const double phi, const double theta, const double psi) {
    phi_ = phi; theta_ = theta; psi_ = psi; }

  /// Set the Euler angles from RotationMatrix.
  void set(const RotationMatrix& matrix);

  /// Return the first angle [-pi, pi] about the z-axis.
  double phi() const { return phi_; }

  /// Return the second angle [0, pi] about the new x-axis.
  double theta() const { return theta_; }

  /// Return the third angle [-pi, pi] about the new z-axis.
  double psi() const { return psi_; }

  /// Compute a RotationMatrix from Euler angles.
  void compute_rotation_matrix(RotationMatrix * matrix) const;

  /// Return true if equal.
  bool is_equal(const Euler& euler, const double tolerance = 1e-8) const;
  
  /// Return human readable format.
  std::string str() const;

 private:
  double phi_;
  double theta_;
  double psi_;
};

}  // namespace feasst

#endif  // FEASST_MATH_EULER_H_
