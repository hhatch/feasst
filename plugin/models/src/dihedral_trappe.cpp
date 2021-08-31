#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "models/include/dihedral_trappe.h"

namespace feasst {

class MapDihedralTRAPPE {
 public:
  MapDihedralTRAPPE() {
    auto obj = MakeDihedralTRAPPE();
    obj->deserialize_map()["DihedralTRAPPE"] = obj;
  }
};

static MapDihedralTRAPPE mapper_ = MapDihedralTRAPPE();

std::shared_ptr<BondFourBody> DihedralTRAPPE::create(std::istream& istr) const {
  return std::make_shared<DihedralTRAPPE>(istr);
}

DihedralTRAPPE::DihedralTRAPPE(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralTRAPPE", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(846 == version, "mismatch version: " << version);
}

void DihedralTRAPPE::serialize_dihedral_trappe_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(846, ostr);
}

void DihedralTRAPPE::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_trappe_(ostr);
}

double DihedralTRAPPE::energy(const double radians, const Bond& dihedral) const {
  FATAL("not implemented");
  return 0.;
}

}  // namespace feasst
