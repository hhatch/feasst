#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "models/include/angle_harmonic.h"

namespace feasst {

class MapAngleHarmonic {
 public:
  MapAngleHarmonic() {
    auto obj = MakeAngleHarmonic();
    obj->deserialize_map()["AngleHarmonic"] = obj;
  }
};

static MapAngleHarmonic mapper_ = MapAngleHarmonic();

std::shared_ptr<BondThreeBody> AngleHarmonic::create(std::istream& istr) const {
  return std::make_shared<AngleHarmonic>(istr);
}

AngleHarmonic::AngleHarmonic(std::istream& istr) : BondThreeBody(istr) {
  // ASSERT(class_name_ == "AngleHarmonic", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4655 == version, "mismatch version: " << version);
}

void AngleHarmonic::serialize_angle_harmonic_(std::ostream& ostr) const {
  serialize_bond_three_body_(ostr);
  feasst_serialize_version(4655, ostr);
}

void AngleHarmonic::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_angle_harmonic_(ostr);
}

double AngleHarmonic::energy(const double radians, const Bond& angle) const {
  //const double equilibrium_degrees = degrees_to_radians(angle.property("equilibrium_degrees"));
  //const double k_energy_per_radian_sq = degrees_to_radians(angle.property("k_energy_per_radian_sq"));
  FATAL("not implemented");
  return 0.;
}

double AngleHarmonic::random_angle_radians(const Angle& angle, const double beta,
    Random * random) const {
  FATAL("not implemented");
  return 0.;
}

}  // namespace feasst
