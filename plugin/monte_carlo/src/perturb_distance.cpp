#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

PerturbDistance::PerturbDistance(argtype args) : PerturbDistance(&args) {
  check_all_used(args);
}
PerturbDistance::PerturbDistance(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbDistance";
  disable_tunable_();
  potential_acceptance_ = integer("potential_acceptance", args, -1);
}

class MapPerturbDistance {
 public:
  MapPerturbDistance() {
    auto obj = MakePerturbDistance();
    obj->deserialize_map()["PerturbDistance"] = obj;
  }
};

static MapPerturbDistance mapper_ = MapPerturbDistance();

std::shared_ptr<Perturb> PerturbDistance::create(std::istream& istr) const {
  return std::make_shared<PerturbDistance>(istr);
}

void PerturbDistance::precompute(TrialSelect * select, System * system) {
  ASSERT(select->has_property("bond_type"), "cannot obtain bond properties");
  bond_type_ = feasst::round(select->property("bond_type"));
}

void PerturbDistance::move(System * system,
                           TrialSelect * select,
                           Random * random) {
  if (potential_acceptance_ == -1) {
    move_once_(system, select, random);
    return;
  }
  int max_attempt = 1e6;
  const double beta = system->thermo_params().beta();
  for (int attempt = 0; attempt < max_attempt; ++attempt) {
    move_once_(system, select, random);
    Potential * poten = system->get_potential(potential_acceptance_);
    const double energy = poten->select_energy(select->mobile(),
                                               system->get_configuration());
    if (random->uniform() < std::exp(-beta*energy)) {
      return;
    }
  }
  FATAL("max_attempt: " << max_attempt << " reached");
}

void PerturbDistance::move_once_(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  Position * site = select->get_mobile()->get_site_position(0, 0);
  DEBUG("mobile " << select->mobile().str());
  DEBUG("old pos " << site->str());
  random->unit_sphere_surface(site);
  double bond_energy;
  site->multiply(random_distance(*system, select, random, &bond_energy));
  select->add_exclude_energy(bond_energy);
  site->add(select->anchor_position(0, 0, *system));
  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDistance::PerturbDistance(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbDistance", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(228 == version, "mismatch version: " << version);
  feasst_deserialize(&bond_type_, istr);
  feasst_deserialize(&potential_acceptance_, istr);
}

void PerturbDistance::serialize_perturb_distance_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(228, ostr);
  feasst_serialize(bond_type_, ostr);
  feasst_serialize(potential_acceptance_, ostr);
}

void PerturbDistance::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
}

double PerturbDistance::random_distance(const System& system,
    const TrialSelect* select,
    Random * random) {
  FATAL("not implemented");
}
double PerturbDistance::random_distance(const System& system,
    const TrialSelect* select,
    Random * random,
    double * bond_energy) {
  const Bond& bond = system.configuration().unique_types().particle(
    select->particle_type()).bond(bond_type_);
  const double beta = system.thermo_params().beta();
  ASSERT(bond_.deserialize_map().count(bond.model()) == 1,
    bond.model() << " not found");
  const BondTwoBody * model = bond_.deserialize_map()[bond.model()].get();
  const double dist = model->random_distance(bond, beta, system.dimension(), random);
  *bond_energy = model->energy(dist, bond);
  return dist;
}

}  // namespace feasst
