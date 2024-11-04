
#ifndef FEASST_CONFINEMENT_ZERO_POTENTIAL_H_
#define FEASST_CONFINEMENT_ZERO_POTENTIAL_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Add the negative of the current energy as Background so that the total energy
  is now zero.
 */
class ZeroPotential : public Action {
 public:
  //@{
  /** @name Arguments
    - configuration_index: index of configuration potential (default: 0).
   */
  explicit ZeroPotential(argtype args = argtype());
  explicit ZeroPotential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<ZeroPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<ZeroPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ZeroPotential(std::istream& istr);
  virtual ~ZeroPotential() {}

  //@}
 private:
  int configuration_index_;
};

inline std::shared_ptr<ZeroPotential> MakeZeroPotential(
    argtype args = argtype()) {
  return std::make_shared<ZeroPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_ZERO_POTENTIAL_H_
