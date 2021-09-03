#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "models/include/dihedral_trappe.h"

namespace feasst {

class MapDihedralTraPPE {
 public:
  MapDihedralTraPPE() {
    auto obj = MakeDihedralTraPPE();
    obj->deserialize_map()["DihedralTraPPE"] = obj;
  }
};

static MapDihedralTraPPE mapper_ = MapDihedralTraPPE();

std::shared_ptr<BondFourBody> DihedralTraPPE::create(std::istream& istr) const {
  return std::make_shared<DihedralTraPPE>(istr);
}

DihedralTraPPE::DihedralTraPPE(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralTraPPE", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(846 == version, "mismatch version: " << version);
}

void DihedralTraPPE::serialize_dihedral_trappe_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(846, ostr);
}

void DihedralTraPPE::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_trappe_(ostr);
}

double DihedralTraPPE::energy(const double radians, const Bond& dihedral) const {
  FATAL("not implemented");
  return 0.;
}

}  // namespace feasst
