
#ifndef FEASST_CHARGE_ELECTRIC_FIELD_H_
#define FEASST_CHARGE_ELECTRIC_FIELD_H_

#include "utils/include/arguments.h"
#include "system/include/model_one_body.h"

namespace feasst {

/**
  Apply an electric field along a dimension.

  \f$U_{ext}(x) = -q E x\f$

  where \f$q\f$ is the charge of a site, \f$E\f$ is the electric field and
  \f$x\f$ is the position along a given dimension.
 */
class ElectricField : public ModelOneBody {
 public:
  //@{
  /** @name Arguments
    - dimension: direction of the electric field (default: 0).
    - field_strength: external electric field strength in force per charge.
      The units are assumed to be Volts/Angstrom, where charge is electron
      charge and energy is in kJ/mol.
   */
  explicit ElectricField(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(const ModelParams& existing) override;

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ElectricField>(istr); }
  explicit ElectricField(std::istream& istr);
  virtual ~ElectricField() {}

  //@}
 private:
  int dimension_;
  double field_strength_;
  double conversion_factor_ = 0.;
};

inline std::shared_ptr<ElectricField> MakeElectricField(
    argtype args = argtype()) {
  return std::make_shared<ElectricField>(args);
}

}  // namespace feasst

#endif  // FEASST_CHARGE_ELECTRIC_FIELD_H_
