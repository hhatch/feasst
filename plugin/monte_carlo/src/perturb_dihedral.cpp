#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/random.h"
#include "math/include/utils_math.h"  // round
#include "monte_carlo/include/perturb_dihedral.h"

namespace feasst {

PerturbDihedral::PerturbDihedral(argtype args)
  : PerturbDihedral(&args) {
  check_all_used(args);
}
PerturbDihedral::PerturbDihedral(argtype * args)
  : PerturbDistanceAngle(args) {
  class_name_ = "PerturbDihedral";
}

class MapPerturbDihedral {
 public:
  MapPerturbDihedral() {
    auto obj = MakePerturbDihedral();
    obj->deserialize_map()["PerturbDihedral"] = obj;
  }
};

static MapPerturbDihedral mapper_ = MapPerturbDihedral();

std::shared_ptr<Perturb> PerturbDihedral::create(std::istream& istr) const {
  return std::make_shared<PerturbDihedral>(istr);
}

void PerturbDihedral::precompute(TrialSelect * select, System * system) {
  PerturbDistanceAngle::precompute(select, system);
  ASSERT(select->has_property("dihedral_type"), "cannot obtain dihedral properties");
  dihedral_type_ = feasst::round(select->property("dihedral_type"));
  DEBUG("dihedral_type_ " << dihedral_type_);
}

double PerturbDihedral::random_dihedral_radians(const System& system,
    const TrialSelect * select, Random * random, double * bond_energy) {
  const Dihedral& dihedral = system.configuration().unique_type(
    select->particle_type()).dihedral(dihedral_type_);
  const double beta = system.thermo_params().beta();
  ASSERT(dihedral_.deserialize_map().count(dihedral.model()) == 1,
    dihedral.model() << " not found");
  const BondFourBody * model = dihedral_.deserialize_map()[dihedral.model()].get();
  const double radians = model->random_dihedral_radians(dihedral, beta, system.dimension(), random);
  *bond_energy += model->energy(radians, dihedral);
  DEBUG("bond_energy " << *bond_energy);
  return radians;
}

void PerturbDihedral::move(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  double bond_energy = 0.;
  const double distance = random_distance(*system, select, random, &bond_energy);
  const double angle = random_angle_radians(*system, select, random, &bond_energy);
  const double dihedral = random_dihedral_radians(*system, select, random, &bond_energy);
  select->add_exclude_energy(bond_energy);
  place_dihedral(distance, angle, dihedral, system, select);
}

void PerturbDihedral::place_dihedral(const double distance,
  const double angle,
  const double dihedral,
  System * system,
  TrialSelect * select) {
  DEBUG("angle " << angle);
  const int dimen = system->configuration().dimension();
  ASSERT(dimen == 3, "not implemented for dimen: " << dimen);
  if (origin_.dimension() == 0) origin_.set_to_origin(dimen);
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  /*
    For given sites j, k and l (anchors), place site, i, according to bond
    length, angle and dihedral.

    l
     \
      k - j
           \
            i

    r_ij = r_i - r_j points from j to i
    - initialize r_ij as r_jk, normalized and multiplied by bond length
    - define the normal n_2 = r_kj cross r_lk = r_jk cross r_kl
    - rotate r_ij by (180 - angle) about n_2
    - rotate r_ij by (180 - dihedral) about r_jk
   */
  const Position& rj = select->anchor_position(0, 0, *system);
  const Position& rk = select->anchor_position(0, 1, *system);
  const Position& rl = select->anchor_position(0, 2, *system);
  DEBUG("rj: " << rj.str() << " rk: " << rk.str() << " rl: " << rl.str());
  rjk_ = rj;
  rjk_.subtract(rk);
  *site = rjk_;

  // normalize site position, centered at rj, and multiply by bond length.
  site->normalize();
  site->multiply(distance);
  rkl_ = rk;
  rkl_.subtract(rl);
  const Position n2 = rjk_.cross_product(rkl_);
  DEBUG("n2 " << n2.str());
  rot_mat_.axis_angle(n2, radians_to_degrees(PI - angle));
  rot_mat_.rotate(origin_, site);
  rot_mat_.axis_angle(rjk_, radians_to_degrees(PI - dihedral));
  rot_mat_.rotate(origin_, site);
  site->add(rj);  // return to frame of reference
  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDihedral::PerturbDihedral(std::istream& istr)
  : PerturbDistanceAngle(istr) {
  ASSERT(class_name_ == "PerturbDihedral", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7579 == version, "mismatch version: " << version);
  feasst_deserialize(&dihedral_type_, istr);
}

void PerturbDihedral::serialize(std::ostream& ostr) const {
  serialize_perturb_dihedral_(ostr);
}

void PerturbDihedral::serialize_perturb_dihedral_(
    std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(7579, ostr);
  feasst_serialize(dihedral_type_, ostr);
}

}  // namespace feasst
