
#ifndef FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_
#define FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(r) = 0 when the angle is between the minimum and maximum specified in
  AngleProperties, otherwise infinity.
 */
class AngleSquareWell : public AngleModel {
 public:
  explicit AngleSquareWell(const argtype& args = argtype()) {}
  double energy(const double theta, const Bond& angle) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AngleSquareWell(std::istream& istr);
  virtual ~AngleSquareWell() {}

 protected:
  void serialize_angle_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<AngleSquareWell> MakeAngleSquareWell(
    const argtype &args = argtype()) {
  return std::make_shared<AngleSquareWell>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_
