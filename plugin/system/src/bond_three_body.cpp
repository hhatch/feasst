
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "system/include/bond_three_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondThreeBody> >& BondThreeBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondThreeBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondThreeBody> >();
  return *ans;
}

void BondThreeBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondThreeBody> BondThreeBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondThreeBody> BondThreeBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondThreeBody::serialize_bond_three_body_(std::ostream& ostr) const {
  feasst_serialize_version(943, ostr);
}

BondThreeBody::BondThreeBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(943 == version, "mismatch version: " << version);
}

double BondThreeBody::energy(const Position& relative01,
    const Position& relative21,
    const Bond& angle) const {
  const double radians = std::acos(relative01.cosine(relative21));
  return energy(radians, angle);
}

void BondThreeBody::random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random) const {
  ASSERT(a2a1m1.model() == a2a1m2.model() &&
         a2a1m1.model() == m1a1m2.model(), "Branch model mismatch: "
    << a2a1m1.model() << " " << a2a1m2.model() << " " << m1a1m2.model());
  const int dimen = 3;
  Position unit_m1(dimen), unit_m2(dimen);
  bool accept = false;
  int attempt = 0;
  while (!accept) {
    random->unit_sphere_surface(&unit_m1);
    random->unit_sphere_surface(&unit_m2);
    const double tm1 = std::acos(unit_m1.coord(0));
    const double tm2 = std::acos(unit_m2.coord(0));
    const double tm12 = std::acos(unit_m1.dot_product(unit_m2));
    if (random->uniform() < std::exp(-beta*(
        energy(tm1, a2a1m1) +
        energy(tm2, a2a1m2) +
        energy(tm12, m1a1m2)))) {
      accept = true;
      *radians_a2a1m1 = tm1;
      *radians_a2a1m2 = tm2;
      *radians_m1a1m2 = tm12;
    }
    ++attempt;
    if (attempt == 1e6) FATAL("max attempts reached");
  }
}

}  // namespace feasst
