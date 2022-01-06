#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt to add a particle.
class TrialTransfer : public TrialFactoryNamed {
 public:
  TrialTransfer(argtype args = argtype());
  TrialTransfer(argtype * args);
//  std::shared_ptr<Trial> create(std::istream& istr) const override {
//    return std::make_shared<TrialTransfer>(istr); }
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransfer>(args); }
//  void serialize(std::ostream& ostr) const override;
//  explicit TrialTransfer(std::istream& istr);
  virtual ~TrialTransfer() {}
};

inline std::shared_ptr<TrialTransfer> MakeTrialTransfer(argtype args = argtype()) {
  return std::make_shared<TrialTransfer>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSFER_H_
