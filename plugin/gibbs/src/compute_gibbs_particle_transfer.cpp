#include <cmath>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select.h"
#include "gibbs/include/compute_gibbs_particle_transfer.h"

namespace feasst {

ComputeGibbsParticleTransfer::ComputeGibbsParticleTransfer() {
  class_name_ = "ComputeGibbsParticleTransfer";
}

class MapComputeGibbsParticleTransfer {
 public:
  MapComputeGibbsParticleTransfer() {
    auto obj = MakeComputeGibbsParticleTransfer();
    obj->deserialize_map()["ComputeGibbsParticleTransfer"] = obj;
  }
};

static MapComputeGibbsParticleTransfer mapper_ = MapComputeGibbsParticleTransfer();

void ComputeGibbsParticleTransfer::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGibbsParticleTransfer");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->add_to_energy_new(criteria->current_energy());
  //acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
  acceptance->add_to_macrostate_shift(1);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  DEBUG("old en " << criteria->current_energy());
  DEBUG("new en " << MAX_PRECISION << acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    acceptance->set_macrostate_shift_type(particle_type);
    const ThermoParams& params = system->thermo_params();
    DEBUG("selprob " << select.probability() << " betamu " << params.beta_mu(particle_type));
    DEBUG("lnselprob " << std::log(select.probability()));
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    acceptance->add_to_ln_metropolis_prob(
      std::log(select.probability())
      + params.beta_mu(particle_type)
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> ComputeGibbsParticleTransfer::create(std::istream& istr) const {
  return std::make_shared<ComputeGibbsParticleTransfer>(istr);
}

ComputeGibbsParticleTransfer::ComputeGibbsParticleTransfer(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeGibbsParticleTransfer", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7345 == version, "mismatch version: " << version);
}

void ComputeGibbsParticleTransfer::serialize_compute_gibbs_particle_transfer_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(7345, ostr);
}

void ComputeGibbsParticleTransfer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_gibbs_particle_transfer_(ostr);
}

}  // namespace feasst
