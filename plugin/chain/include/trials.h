#ifndef FEASST_CHAIN_TRIALS_H_
#define FEASST_CHAIN_TRIALS_H_

#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Swap the types of two sites in a particle.
  args:
  - site_type1: type of site to swap.
  - site_type2: type of other site to swap.
 */
std::shared_ptr<Trial> MakeTrialSwapSites(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIALS_H_

