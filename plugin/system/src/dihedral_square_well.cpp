#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "system/include/dihedral_square_well.h"

namespace feasst {

class MapDihedralSquareWell {
 public:
  MapDihedralSquareWell() {
    auto obj = MakeDihedralSquareWell();
    obj->deserialize_map()["DihedralSquareWell"] = obj;
  }
};

static MapDihedralSquareWell mapper_ = MapDihedralSquareWell();

std::shared_ptr<BondFourBody> DihedralSquareWell::create(std::istream& istr) const {
  return std::make_shared<DihedralSquareWell>(istr);
}

DihedralSquareWell::DihedralSquareWell(std::istream& istr) : BondFourBody(istr) {
  // ASSERT(class_name_ == "DihedralSquareWell", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2947 == version, "mismatch version: " << version);
}

void DihedralSquareWell::serialize_dihedral_square_well_(std::ostream& ostr) const {
  serialize_bond_four_body_(ostr);
  feasst_serialize_version(2947, ostr);
}

void DihedralSquareWell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_dihedral_square_well_(ostr);
}

double DihedralSquareWell::energy(
    const Position& ri,
    const Position& rj,
    const Position& rk,
    const Position& rl,
    const Dihedral& dihedral) const {
  return 0.;
}

}  // namespace feasst
