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
  const int dihedral_type = feasst::round(select->property("dihedral_type"));
  const Dihedral& dihedral = system->configuration().unique_type(
    select->particle_type()).dihedral(dihedral_type);
  dihedral_ = degrees_to_radians(dihedral.property("theta0"));
  DEBUG("dihedral_ " << dihedral_);
  const int dimen = system->configuration().dimension();
  ASSERT(dimen == 3, "not implemented for dimen " << dimen);
  if (dihedral.has_property("spring_constant")) {
    spring_constant_ = dihedral.property("spring_constant");
  }
}

double PerturbDihedral::random_dihedral(Random * random,
    const double beta,
    const int dimension) const {
  if (is_rigid()) return dihedral_;
  FATAL("for trappe, different functional form to select based on U");
  return random->bond_angle(dihedral_, beta*spring_constant_, 2, dimension);
}

void PerturbDihedral::move(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  const double beta = system->thermo_params().beta();
  const double distance = random_distance(*system, select, random);
  const double angle = random_angle(random, beta, system->dimension());
  const double dihedral = random_dihedral(random, beta, system->dimension());
  place_dihedral(distance, angle, dihedral, system, select);
}

void PerturbDihedral::place_dihedral(const double distance,
  const double angle,
  const double dihedral,
  System * system,
  TrialSelect * select) {
  const int dimen = system->configuration().dimension();
  ASSERT(dimen == 3, "not implemented for dimen: " << dimen);
  if (origin_.dimension() == 0) origin_.set_to_origin(dimen);
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  /*
    For given sites j, k (anchors), place site, i, according to bond angle,
    length and dihedral.

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
  const Position& rl = select->anchor_position(0, 1, *system);
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
  feasst_deserialize(&dihedral_, istr);
  feasst_deserialize(&spring_constant_, istr);
}

void PerturbDihedral::serialize(std::ostream& ostr) const {
  serialize_perturb_dihedral_(ostr);
}

void PerturbDihedral::serialize_perturb_dihedral_(
    std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(7579, ostr);
  feasst_serialize(dihedral_, ostr);
  feasst_serialize(spring_constant_, ostr);
}

bool PerturbDihedral::is_rigid() const {
  if (std::abs(spring_constant_ + 1) < NEAR_ZERO) return true;
  return false;
}

}  // namespace feasst
