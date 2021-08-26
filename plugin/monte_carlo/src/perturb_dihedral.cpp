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
  return random->bond_angle(dihedral_, beta*spring_constant_, 2, dimension);
}

void PerturbDihedral::move(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  const double distance = random_distance(random,
    system->thermo_params().beta(),
    system->dimension());
  const double angle = random_angle(random, system->thermo_params().beta(),
    system->dimension());
  DEBUG("angle: " << angle);
  place_in_circle(distance, angle, system, select, random);
}

void PerturbDihedral::place_in_circle(const double distance,
  const double angle,
  System * system,
  TrialSelect * select,
  Random * random) {
  const int dimension = system->configuration().dimension();
  if (origin_.dimension() == 0) origin_.set_to_origin(dimension);
  if (orthogonal_jk_.dimension() == 0) orthogonal_jk_.set_to_origin(dimension);
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  /*
    For given sites j, k (anchors), place site, i, according to bond angle and
    length

    k - j
         \
   angle  i

    r_jk = r_j - r_k points from k to j
    - set the bond length by normalizing r_jk and multiplying by bond length
    - set the bond angle by rotating r_jk by 180-angle about axis orthogonal
      to r_jk.
   */
  // set site to the vector r_jk = r_j - r_k and store this vector
  const Position& rj = select->anchor_position(0, 0, *system);
  DEBUG("rj: " << rj.str());
  const Position& rk = select->anchor_position(0, 1, *system);
  DEBUG("rk: " << rk.str());
  rjk_ = rj;
  rjk_.subtract(rk);
  *site = rjk_;
  DEBUG("rjk: " << rjk_.str());

  // normalize site position, centered at rj, and multiply by bond length.
  site->normalize();
  site->multiply(distance);
  DEBUG("rjk norm*L: " << rjk_.str());

  // rotate site by (PI-bond_angle). If 3D, about vector orthogonal to r_jk.
  if (dimension == 3) orthogonal_jk_.orthogonal(*site);
  DEBUG("ortho " << orthogonal_jk_.str());
  rot_mat_.axis_angle(orthogonal_jk_, radians_to_degrees(PI - angle));
  DEBUG("site == rjk: " << site->str());
  rot_mat_.rotate(origin_, site);
  DEBUG("site rotated to angle: " << site->str());

  // If 3D, randomly spin site about rjk.
  if (dimension == 3) {
    rot_mat_.axis_angle(rjk_, 360.*random->uniform());
    rot_mat_.rotate(origin_, site);
  }

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
