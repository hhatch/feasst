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
  INFO("ComputeGibbsParticleTransfer");
  INFO("lnmet " << acceptance->ln_metropolis_prob());

  // del
  std::vector<TrialStage*> del_stages = {(*stages)[1]};
  compute_rosenbluth(1, criteria, system, acceptance, &del_stages, random);
  int config_del = 0;
  int config_add = 1;
  if ((*stages)[0]->trial_select().configuration_index() == 0) {
    config_add = 0;
    config_del = 1;
  }
  INFO("config_add " << config_add);
  INFO("config_del " << config_del);
  acceptance->set_energy_new(criteria->current_energy(config_del) - acceptance->energy_old(config_del), config_del);
  acceptance->set_energy_profile_new(criteria->current_energy_profile(config_del), config_del);
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config_del), config_del);
  acceptance->add_to_macrostate_shift(-1, config_del);

  INFO("energy contribution of config " << config_del << " particle to be deleted: " << acceptance->energy_old(config_del));

  INFO("en 0 old " << criteria->current_energy(0));
  INFO("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
  INFO("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));

  // add
  std::vector<TrialStage*> add_stages = {(*stages)[0]};
  compute_rosenbluth(0, criteria, system, acceptance, &add_stages, random);
  acceptance->add_to_energy_new(criteria->current_energy(config_add), config_add);
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config_add), config_add);
  acceptance->add_to_macrostate_shift(1, config_add);

  INFO("energy contribution of config " << config_add << " particle to be added: " << acceptance->energy_new(config_add));

  INFO("lnmet " << acceptance->ln_metropolis_prob());
//  INFO("en 0 old " << criteria->current_energy(0));
//  INFO("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
//  INFO("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));
//  INFO("en 1 old " << criteria->current_energy(1));
//  INFO("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
//  INFO("en 1 old acc " << MAX_PRECISION << acceptance->energy_old(1));
  { // Metropolis
//    const Configuration& config = system->configuration();
//    const TrialSelect& select = (*stages)[0]->trial_select();
//    const int particle_index = select.mobile().particle_index(0);
//    const int particle_type = config.select_particle(particle_index).type();
//    acceptance->set_macrostate_shift_type(particle_type, config_add);
//    acceptance->set_macrostate_shift_type(particle_type, config_del);
//    INFO("lnselprob " << std::log(select.probability()));
//    INFO("lnmet " << acceptance->ln_metropolis_prob());
//    acceptance->add_to_ln_metropolis_prob(
//      std::log(select.probability())
//    );
    INFO("lnmet " << acceptance->ln_metropolis_prob());
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
