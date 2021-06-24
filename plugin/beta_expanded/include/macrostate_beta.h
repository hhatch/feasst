
#ifndef FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
#define FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the inverse temperature, \f$\beta\f$.
 */
class MacrostateBeta : public Macrostate {
 public:
  // HWH consider depreciating this interface
  MacrostateBeta(const Histogram& histogram, argtype args = argtype());
  MacrostateBeta(const Histogram& histogram, argtype * args);
  explicit MacrostateBeta(argtype args = argtype());
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override {
    return system.thermo_params().beta(); }
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateBeta(std::istream& istr);
  virtual ~MacrostateBeta() {}

 private:
  const std::string class_name_ = "MacrostateBeta";
};

inline std::shared_ptr<MacrostateBeta> MakeMacrostateBeta(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostateBeta>(histogram, args);
}

inline std::shared_ptr<MacrostateBeta> MakeMacrostateBeta(
    argtype args = argtype()) {
  return std::make_shared<MacrostateBeta>(args);
}

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
