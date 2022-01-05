#ifndef FEASST_MONTE_CARLO_TRIALS_H_
#define FEASST_MONTE_CARLO_TRIALS_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Attempt TrialAdd or TrialRemove with equal probability.
std::shared_ptr<TrialFactory> MakeTrialTransfer(
  argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIALS_H_
