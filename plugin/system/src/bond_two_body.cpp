
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "system/include/bond_two_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondTwoBody> >& BondTwoBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondTwoBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondTwoBody> >();
  return *ans;
}

void BondTwoBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondTwoBody> BondTwoBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondTwoBody> BondTwoBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondTwoBody::serialize_bond_two_body_(std::ostream& ostr) const {
  feasst_serialize_version(264, ostr);
}

BondTwoBody::BondTwoBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(264 == version, "mismatch version: " << version);
}

double BondTwoBody::energy(const Position& relative, const Bond& bond) const {
  return energy(relative.distance(), bond);
}

//double BondTwoBody::random_distance(const Bond& bond, const double beta,
//    Random * random) const {
//  double length = random->uniform_real(0, 2*bond.property("equilibrium_length");
//  int attempt = 0;
//  while (attempt < 1e6) {
//  }
//  FATAL("max attempts reached");
//}

}  // namespace feasst
