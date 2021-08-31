
#ifndef FEASST_MODELS_BOND_HARMONIC_H_
#define FEASST_MODELS_BOND_HARMONIC_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  U(r) = k_energy_per_length_sq*(r - equilibrium_length)^2
  with parameters given in Bond Properties

  Note that the optimized Gaussian implementation may assume 3D.
 */
class BondHarmonic : public BondLength {
 public:
  BondHarmonic() {}
  double energy(const double distance, const Bond& bond) const override;
  double random_distance(const Bond& bond, const double beta,
    Random * random) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit BondHarmonic(std::istream& istr);
  virtual ~BondHarmonic() {}

 protected:
  void serialize_bond_harmonic_(std::ostream& ostr) const;
};

inline std::shared_ptr<BondHarmonic> MakeBondHarmonic() {
  return std::make_shared<BondHarmonic>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_BOND_HARMONIC_H_
