#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to add a particle.
class TrialRemove : public Trial {
 public:
  TrialRemove(argtype args = argtype());
  TrialRemove(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemove>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemove>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemove(std::istream& istr);
  virtual ~TrialRemove() {}
};

inline std::shared_ptr<TrialRemove> MakeTrialRemove(argtype args = argtype()) {
  return std::make_shared<TrialRemove>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_REMOVE_H_
