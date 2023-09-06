#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"

namespace feasst {

ComputeGibbsVolumeTransfer::ComputeGibbsVolumeTransfer() {
  class_name_ = "ComputeGibbsVolumeTransfer";
}

class MapComputeGibbsVolumeTransfer {
 public:
  MapComputeGibbsVolumeTransfer() {
    auto obj = MakeComputeGibbsVolumeTransfer();
    obj->deserialize_map()["ComputeGibbsVolumeTransfer"] = obj;
  }
};

static MapComputeGibbsVolumeTransfer mapper_ = MapComputeGibbsVolumeTransfer();

void ComputeGibbsVolumeTransfer::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGibbsVolumeTransfer");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());

  // del
  std::vector<TrialStage*> del_stages = {(*stages)[1]};
  compute_rosenbluth(1, criteria, system, acceptance, &del_stages, random);
  int config_del = 0;
  int config_add = 1;
  if ((*stages)[0]->trial_select().configuration_index() == 0) {
    config_add = 0;
    config_del = 1;
  }
  DEBUG("config_add " << config_add);
  DEBUG("config_del " << config_del);
  acceptance->set_energy_new(criteria->current_energy(config_del) - acceptance->energy_old(config_del), config_del);
  acceptance->set_energy_profile_new(criteria->current_energy_profile(config_del), config_del);
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config_del), config_del);
  acceptance->add_to_macrostate_shift(-1, config_del);

  DEBUG("energy contribution of config " << config_del << " particle to be deleted: " << acceptance->energy_old(config_del));

  DEBUG("en 0 old " << criteria->current_energy(0));
  DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
  DEBUG("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));

  // add
  std::vector<TrialStage*> add_stages = {(*stages)[0]};
  compute_rosenbluth(0, criteria, system, acceptance, &add_stages, random);
  acceptance->add_to_energy_new(criteria->current_energy(config_add), config_add);
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config_add), config_add);
  acceptance->add_to_macrostate_shift(1, config_add);

  DEBUG("energy contribution of config " << config_add << " particle to be added: " << acceptance->energy_new(config_add));

  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
//  DEBUG("en 0 old " << criteria->current_energy(0));
//  DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
//  DEBUG("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));
//  DEBUG("en 1 old " << criteria->current_energy(1));
//  DEBUG("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
//  DEBUG("en 1 old acc " << MAX_PRECISION << acceptance->energy_old(1));
  { // Metropolis
    const Configuration& conf_add = system->configuration(config_add);
    const Configuration& conf_del = system->configuration(config_del);
    const TrialSelect& select_add = (*stages)[0]->select();
    //const TrialSelect& select_del = (*stages)[1]->select();
    const int particle_add = select_add.mobile().particle_index(0);
    //const int particle_del = select_del.mobile().particle_index(0);
    const int particle_type = conf_add.select_particle(particle_add).type();
    acceptance->set_macrostate_shift_type(particle_type, config_add);
    acceptance->set_macrostate_shift_type(particle_type, config_del);
//    DEBUG("lnselprob " << std::log(select.probability()));
//    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    const int num_particles_from_add = conf_add.num_particles_of_type(particle_type);
    const int num_particles_from_del = conf_del.num_particles_of_type(particle_type);
    const double vol_from_add = conf_add.domain().volume();
    const double vol_from_del = conf_del.domain().volume();
    acceptance->add_to_ln_metropolis_prob(
      num_particles_from_del/static_cast<double>(num_particles_from_add+1)*
      vol_from_add/vol_from_del
      //std::log(select.probability())
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> ComputeGibbsVolumeTransfer::create(std::istream& istr) const {
  return std::make_shared<ComputeGibbsVolumeTransfer>(istr);
}

ComputeGibbsVolumeTransfer::ComputeGibbsVolumeTransfer(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeGibbsVolumeTransfer", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7345 == version, "mismatch version: " << version);
}

void ComputeGibbsVolumeTransfer::serialize_compute_gibbs_volume_transfer_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(7345, ostr);
}

void ComputeGibbsVolumeTransfer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_gibbs_volume_transfer_(ostr);
}

}  // namespace feasst
