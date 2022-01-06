#ifndef FEASST_CHAIN_TRIALS_H_
#define FEASST_CHAIN_TRIALS_H_

#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/// Rigidly pivot an end segment of a chain to a random orientation.
std::shared_ptr<Trial> MakeTrialPivot(argtype args = argtype());

/**
  Reptate a linear chain by taking one end and adding it to the other end.
  For heteropolymers, this perturbation changes the composition of all
  particles with the same type.
  Thus, individual heteropolymers should be added as unique particle types.
  The bond length is taken as the bond between site 0 and 1, and is assumed
  to be constant.
  Thus, as currently implemented, heteropolymers must have the same bond length.
  This trial may not be compatible with angle and dihedral potentials.
  Instead, use "reptate" in TrialGrow.
 */
std::shared_ptr<Trial> MakeTrialReptate(argtype args = argtype());

/**
  Swap the types of two sites in a particle.
  args:
  - site_type1: type of site to swap.
  - site_type2: type of other site to swap.
 */
std::shared_ptr<Trial> MakeTrialSwapSites(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIALS_H_

