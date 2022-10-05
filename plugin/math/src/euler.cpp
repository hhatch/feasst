#include <cmath>
#include <vector>
#include <sstream>
#include "utils/include/serialize.h"
#include "math/include/euler.h"

namespace feasst {

void Euler::compute_rotation_matrix(RotationMatrix * matrix) const {
  if (matrix->num_rows() == 0) {
    matrix->set_size(3, 3);
  }
  const double sphi = std::sin(phi_);
  const double cphi = std::cos(phi_);
  const double stheta = std::sin(theta_);
  const double ctheta = std::cos(theta_);
  const double spsi = std::sin(psi_);
  const double cpsi = std::cos(psi_);
  // See https://mathworld.wolfram.com/EulerAngles.html.
  matrix->set_value(0, 0, cpsi*cphi - ctheta*sphi*spsi); //a11
  matrix->set_value(0, 1, cpsi*sphi + ctheta*cphi*spsi); //a12
  matrix->set_value(0, 2, spsi*stheta);                  //a13
  matrix->set_value(1, 0, -spsi*cphi - ctheta*sphi*cpsi);//a21
  matrix->set_value(1, 1, -spsi*sphi + ctheta*cphi*cpsi);//a22
  matrix->set_value(1, 2, cpsi*stheta);                  //a23
  matrix->set_value(2, 0, stheta*sphi);                  //a31
  matrix->set_value(2, 1, -stheta*cphi);                 //a32
  matrix->set_value(2, 2, ctheta);                       //a33
}

void Euler::set(const RotationMatrix& matrix) {
  const std::vector<std::vector<double> >& mat = matrix.matrix();
  if (std::abs(mat[2][2] - 1.) < 1e-12) {
    theta_ = std::acos(1.);
  } else {
    theta_ = std::acos(mat[2][2]);
  }
  const double stheta = std::sin(theta_);
  if (std::abs(stheta) < 1e-8) {
    phi_ = 0.;
    psi_ = std::atan2(mat[0][1], mat[0][0]);
  } else {
    phi_ = std::atan2(mat[2][0]/stheta,-mat[2][1]/stheta);
    psi_ = std::atan2(mat[0][2]/stheta, mat[1][2]/stheta);
  }
}

bool Euler::is_equal(const Euler& euler, const double tolerance) const {
  if (std::abs(phi_ - euler.phi()) > tolerance) {
    return false;
  } else if (std::abs(theta_ - euler.theta()) > tolerance) {
    return false;
  } else if (std::abs(psi_ - euler.psi()) > tolerance) {
    return false;
  }
  return true;
}

std::string Euler::str() const {
  std::stringstream ss;
  ss << phi_ << "," << theta_ << "," << psi_;
  return ss.str();
}

void Euler::serialize(std::ostream& ostr) const {
  feasst_serialize_version(3509, ostr);
  feasst_serialize(phi_, ostr);
  feasst_serialize(theta_, ostr);
  feasst_serialize(psi_, ostr);
}

Euler::Euler(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3509, "unrecognized version: " << version);
  feasst_deserialize(&phi_, istr);
  feasst_deserialize(&theta_, istr);
  feasst_deserialize(&psi_, istr);
}

}  // namespace feasst
