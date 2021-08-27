
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "system/include/bond_four_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondFourBody> >& BondFourBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondFourBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondFourBody> >();
  return *ans;
}

void BondFourBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondFourBody> BondFourBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondFourBody> BondFourBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondFourBody::serialize_bond_four_body_(std::ostream& ostr) const {
  feasst_serialize_version(7509, ostr);
}

BondFourBody::BondFourBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(7509 == version, "mismatch version: " << version);
}

}  // namespace feasst
