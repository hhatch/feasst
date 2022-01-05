#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialRotate(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbRotate>(&args),
    "TrialRotate",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransfer(argtype args) {
  argtype orig_args = args;
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAdd(orig_args));
  factory->add(MakeTrialRemove(orig_args));
  return factory;
}

std::shared_ptr<Trial> MakeTrialVolume(argtype args) {
  auto trial = MakeTrial(&args);
  trial->set_description("TrialVolume");
  trial->add_stage(
    std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbVolume>(&args),
    &args);
  trial->set(MakeTrialComputeVolume());
  return trial;
}

}  // namespace feasst
