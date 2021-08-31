
#ifndef FEASST_SYSTEM_DIHEDRAL_SQUARE_WELL_H_
#define FEASST_SYSTEM_DIHEDRAL_SQUARE_WELL_H_

#include <memory>
#include "system/include/bond_four_body.h"

namespace feasst {

/**
 */
class DihedralSquareWell : public BondFourBody {
 public:
  explicit DihedralSquareWell(const argtype& args = argtype()) {}
  double energy(
      const Position& ri,
      const Position& rj,
      const Position& rk,
      const Position& rl,
      const Dihedral& dihedral) const override;
  std::shared_ptr<BondFourBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit DihedralSquareWell(std::istream& istr);
  virtual ~DihedralSquareWell() {}

 protected:
  void serialize_dihedral_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<DihedralSquareWell> MakeDihedralSquareWell(
    const argtype &args = argtype()) {
  return std::make_shared<DihedralSquareWell>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_DIHEDRAL_SQUARE_WELL_H_
