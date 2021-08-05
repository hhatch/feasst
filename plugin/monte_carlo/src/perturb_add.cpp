#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_add.h"

namespace feasst {

PerturbAdd::PerturbAdd(argtype args) : PerturbAdd(&args) {
  check_all_used(args);
}
PerturbAdd::PerturbAdd(argtype * args) : Perturb(args) {
  class_name_ = "PerturbAdd";
  delay_add_ = boolean("delay_add", args, true);
  disable_tunable_();
}

class MapPerturbAdd {
 public:
  MapPerturbAdd() {
    auto obj = MakePerturbAdd();
    obj->deserialize_map()["PerturbAdd"] = obj;
  }
};

static MapPerturbAdd mapper_ = MapPerturbAdd();

std::shared_ptr<Perturb> PerturbAdd::create(std::istream& istr) const {
  return std::make_shared<PerturbAdd>(istr);
}

void PerturbAdd::add(
    System * system,
    TrialSelect * select,
    Random * random,
    const Position& center,
    const bool is_position_held) {
  DEBUG("is_position_held " << is_position_held);
  DEBUG(select->mobile().str());
  Configuration* config = system->get_configuration();
  if (!delay_add_) config->revive(select->mobile());
  const int particle_type = config->select_particle(
    select->mobile().particle_index(0)
  ).type();
  DEBUG("type " << particle_type);
  for (const Select& ghost : config->ghosts()) {
    DEBUG("ghost " << ghost.str());
  }
  if (center.dimension() == 0) {
    anywhere_.perturb(system, select, random, is_position_held);
  } else {
    anywhere_.set_position(center, system, select);
  }
  set_revert_possible(true, select);
  set_finalize_possible(true, select);

  // setting trial state should go last so other perturbs do not overwrite
  DEBUG("setting trial state 3");
  select->set_trial_state(3);
}

void PerturbAdd::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
//    DEBUG(revert_select()->mobile().str());
//    DEBUG("nump " << system->configuration().num_particles());
//    system->revert(revert_select()->mobile());
    anywhere_.revert(system);
    if (!delay_add_) {
      system->get_configuration()->remove_particles(revert_select()->mobile());
    }
  }
}

void PerturbAdd::finalize(System * system) {
  DEBUG("finalize_possible " << finalize_possible());
  if (finalize_possible()) {
//    INFO("finalizing mobile " << finalize_select()->mobile().str());
//    system->finalize(finalize_select()->mobile());
    if (delay_add_) {
      system->get_configuration()->revive(finalize_select()->mobile());
    }
  }
}

std::string PerturbAdd::status_header() const {
  std::stringstream ss;
  return ss.str();
}

std::string PerturbAdd::status() const {
  std::stringstream ss;
  return ss.str();
}

PerturbAdd::PerturbAdd(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(730 == version, "mismatch version: " << version);
  feasst_deserialize(&delay_add_, istr);
}

void PerturbAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(730, ostr);
  feasst_serialize(delay_add_, ostr);
}

}  // namespace feasst
