
#ifndef FEASST_SYSTEM_ANGLE_HARMONIC_H_
#define FEASST_SYSTEM_ANGLE_HARMONIC_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(angle) = k_energy_per_radian_sq*(angle - equilibrium_degrees)^2
  with parameters given in Angle Properties

  Note that the optimized Gaussian implementation may assume 3D.
 */
class AngleHarmonic : public BondThreeBody {
 public:
  AngleHarmonic() {}
  double energy(const double radians, const Bond& angle) const override;
  double random_angle_radians(const Angle& angle, const double beta,
    Random * random) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AngleHarmonic(std::istream& istr);
  virtual ~AngleHarmonic() {}

 protected:
  void serialize_angle_harmonic_(std::ostream& ostr) const;
};

inline std::shared_ptr<AngleHarmonic> MakeAngleHarmonic() {
  return std::make_shared<AngleHarmonic>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ANGLE_HARMONIC_H_
