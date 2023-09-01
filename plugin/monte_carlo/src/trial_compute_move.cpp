#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

TrialComputeMove::TrialComputeMove(argtype args) : TrialComputeMove(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TrialComputeMove::TrialComputeMove(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeMove";
}

class MapTrialComputeMove {
 public:
  MapTrialComputeMove() {
    auto obj = MakeTrialComputeMove();
    obj->deserialize_map()["TrialComputeMove"] = obj;
  }
};

static MapTrialComputeMove mapper_ = MapTrialComputeMove();

void TrialComputeMove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  INFO("TrialComputeMove");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  INFO("New");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const int config = stages->front()->select().configuration_index();
  INFO("config " << config);
  acceptance->set_configuration_index(config);
  INFO("current en: " << criteria->current_energy(config));
  INFO("old en: " << acceptance->energy_old());
  INFO("new en: " << acceptance->energy_new());
  INFO("energy change: " << acceptance->energy_new() - acceptance->energy_old());
  INFO("acceptance energy profile old " << feasst_str(acceptance->energy_profile_old()));
  if ((*stages)[0]->is_new_only()) {
    //acceptance->set_energy_new(acceptance->energy_new());
  } else {
    const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
    acceptance->set_energy_new(criteria->current_energy(config) + delta_energy);
    INFO("new profile bf " << feasst_str(acceptance->energy_profile_new()));
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config));
    INFO("new profile af " << feasst_str(acceptance->energy_profile_new()));
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());
    INFO("new profile af2 " << feasst_str(acceptance->energy_profile_new()));
  }
}

std::shared_ptr<TrialCompute> TrialComputeMove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeMove>(istr);
}

TrialComputeMove::TrialComputeMove(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void TrialComputeMove::serialize_trial_compute_move_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void TrialComputeMove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_move_(ostr);
}

}  // namespace feasst
