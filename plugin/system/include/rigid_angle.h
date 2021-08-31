
#ifndef FEASST_SYSTEM_RIGID_ANGLE_H_
#define FEASST_SYSTEM_RIGID_ANGLE_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(r) = 0 when the absolute value of the difference between the actual angle
  and the degrees parameter is less than delta.
  Otherwise, U(r) = NEAR_INFINITY.
 */
class RigidAngle : public AngleModel {
 public:
  explicit RigidAngle(const argtype& args = argtype()) {}
  double energy(const double theta, const Bond& angle) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit RigidAngle(std::istream& istr);
  virtual ~RigidAngle() {}

 protected:
  void serialize_rigid_angle_(std::ostream& ostr) const;
};

inline std::shared_ptr<RigidAngle> MakeRigidAngle(
    const argtype &args = argtype()) {
  return std::make_shared<RigidAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_RIGID_ANGLE_H_
