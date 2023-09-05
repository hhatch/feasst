#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_compute_remove.h"

namespace feasst {

TrialComputeRemove::TrialComputeRemove(argtype args) : TrialComputeRemove(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TrialComputeRemove::TrialComputeRemove(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeRemove";
}

class MapTrialComputeRemove {
 public:
  MapTrialComputeRemove() {
    auto obj = MakeTrialComputeRemove();
    obj->deserialize_map()["TrialComputeRemove"] = obj;
  }
};

static MapTrialComputeRemove mapper_ = MapTrialComputeRemove();

void TrialComputeRemove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeRemove");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  const int iconf = stages->front()->select().configuration_index();
  acceptance->set_energy_new(criteria->current_energy(iconf) - acceptance->energy_old(iconf), iconf);
  acceptance->set_energy_profile_new(criteria->current_energy_profile(iconf), iconf);
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(iconf), iconf);
  acceptance->add_to_macrostate_shift(-1, iconf);
  { // Metropolis
    const Configuration& config = system->configuration(iconf);
    const double volume = config.domain().volume();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    acceptance->set_macrostate_shift_type(particle_type, iconf);
    DEBUG("volume " << volume << " selprob " << select.probability() << " betamu " << system->thermo_params().beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      - std::log(volume*select.probability())
      - system->thermo_params().beta_mu(particle_type)
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> TrialComputeRemove::create(
    std::istream& istr) const {
  return std::make_shared<TrialComputeRemove>(istr);
}

TrialComputeRemove::TrialComputeRemove(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(332 == version, "mismatch version: " << version);
}

void TrialComputeRemove::serialize_trial_compute_remove_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(332, ostr);
}

void TrialComputeRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_remove_(ostr);
}

}  // namespace feasst
