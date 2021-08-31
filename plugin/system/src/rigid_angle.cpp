#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "system/include/rigid_angle.h"

namespace feasst {

class MapRigidAngle {
 public:
  MapRigidAngle() {
    auto obj = MakeRigidAngle();
    obj->deserialize_map()["RigidAngle"] = obj;
  }
};

static MapRigidAngle mapper_ = MapRigidAngle();

std::shared_ptr<BondThreeBody> RigidAngle::create(std::istream& istr) const {
  return std::make_shared<RigidAngle>(istr);
}

RigidAngle::RigidAngle(std::istream& istr) : AngleModel(istr) {
  // ASSERT(class_name_ == "RigidAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5968 == version, "mismatch version: " << version);
}

void RigidAngle::serialize_rigid_angle_(std::ostream& ostr) const {
  serialize_bond_three_body_(ostr);
  feasst_serialize_version(5968, ostr);
}

void RigidAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_rigid_angle_(ostr);
}

double RigidAngle::energy(const double theta, const Bond& angle) const {
  const double radians = degrees_to_radians(angle.property("degrees"));
  const double delta = degrees_to_radians(angle.property("delta"));
  ASSERT(!std::isnan(radians), "radians is nan");
  if (std::abs(theta - radians) > delta) {
    return NEAR_INFINITY;
  }
  return 0.;
}

}  // namespace feasst
