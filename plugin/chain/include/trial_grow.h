#ifndef FEASST_CHAIN_TRIAL_GROW_H_
#define FEASST_CHAIN_TRIAL_GROW_H_

#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Manually describe partial or full-particle growth using configurational bias.
  The input is a vector of argtype, where each argtype represents a TrialStage
  for growing one (or rarely, multiple) sites.

  The following options may only be used in the first argtype.
  - particle_type: type of particle in Configuration (always required).
  - site: site index in particle_type to transfer/regrow (always required).
  - weight: weight of selection of this trial (default: see Trial).
  - transfer: if true, create add and remove trial with equal weight
    (default: false).
  - regrow: if true, place anywhere in the domain (default: false).
  - transfer_avb: if true, same as transfer but with AVB Add/Remove for the
    first stage (default: false).
  - regrow_avb2: if true, regrow using AVB2 in the first stage (default: false).
  - regrow_avb4: if true, regrow using AVB4 in the first stage (default: false).
  - translate: if true (default: false), translate site (which is required arg
    for TrialSelectParticle).
    In addition, must have number of stages equal to number of sites.

  The following options may be used in any argtype.
  If used in the first, then its a partial regrowth move.
  - bond: if used, add TrialSelectBond and PerturbDistance.
    Requires arguments described in TrialSelectBond.
  - angle: if used, adds TrialSelectAngle and PerturbDistanceAngle.
    Requires arguments described in TrialSelectAngle and TrialSelectBond.
  - dihedral: if used, adds TrialSelectDihedral and PerturbDihedral.
    Requires arguments described in SelectDihedral, Angle and Bond.
  - branch: if used, adds SelectBranch and PerturbBranch.
    Requires arguments described in SelectBranch, Angle and Bond.
  - TrialStage arguments: num_steps, reference_index, new_only, etc.

  Note that only one of bond, angle or branch may be true for a given stage.
 */
std::shared_ptr<TrialFactory> MakeTrialGrow(std::vector<argtype> args,
  /// Optionally, set the default values for the following TrialStage arguments:
  /// num_steps, reference_index and new_only.
  /// Any option applied by the above args overwrites this option.
  const argtype& default_args = argtype());

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_GROW_H_
