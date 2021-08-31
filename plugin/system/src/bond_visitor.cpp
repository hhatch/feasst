#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/configuration.h"
#include "system/include/bond_visitor.h"

namespace feasst {

class MapBondVisitor {
 public:
  MapBondVisitor() {
    auto obj = MakeBondVisitor();
    obj->deserialize_map()["BondVisitor"] = obj;
  }
};

static MapBondVisitor mapper_ = MapBondVisitor();

std::map<std::string, std::shared_ptr<BondVisitor> >& BondVisitor::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondVisitor> >* ans =
     new std::map<std::string, std::shared_ptr<BondVisitor> >();
  return *ans;
}

BondVisitor::BondVisitor(argtype args) {
  verbose_ = boolean("verbose", &args, false);
  if (VERBOSE_LEVEL == 5) verbose_ = true;
  check_all_used(args);
}

void BondVisitor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_visitor_(ostr);
}

std::shared_ptr<BondVisitor> BondVisitor::create(std::istream& istr) const {
  return std::make_shared<BondVisitor>(istr);
}

std::shared_ptr<BondVisitor> BondVisitor::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondVisitor::serialize_bond_visitor_(std::ostream& ostr) const {
  feasst_serialize_version(303, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize(verbose_, ostr);
}

BondVisitor::BondVisitor(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(303 == version, "mismatch version: " << version);
  feasst_deserialize(&energy_, istr);
  feasst_deserialize(&verbose_, istr);
}

void BondVisitor::compute_two(const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Bond& bond : config.particle_type(part_type).bonds()) {
      const Position& position0 = part.site(bond.site(0)).position();
      const Position& position1 = part.site(bond.site(1)).position();
      Position relative = position0;
      relative.subtract(position1);
      const Bond& bond_type = config.unique_type(part_type).bond(bond.type());
      ASSERT(bond_.deserialize_map().count(bond_type.model()) == 1,
        "bond model " << bond_type.model() << " not recognized.");
      en += bond_.deserialize_map()[bond_type.model()]->energy(
        relative, bond_type);
      if (verbose_) {
        if (std::abs(en) > NEAR_ZERO) {
          INFO("bond ij " << part_index << " " << bond.site(0) << " "
            << bond.site(1) << " sq " << relative.squared_distance());
        }
      }
    }
  }
  energy_two_body_ = en;
}

void BondVisitor::compute_three(
    const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Angle& angle : config.particle_type(part_type).angles()) {
      const Position& position0 = part.site(angle.site(0)).position();
      const Position& position1 = part.site(angle.site(1)).position();
      const Position& position2 = part.site(angle.site(2)).position();
      Position relative01 = position0;
      Position relative21 = position2;
      relative01.subtract(position1);
      relative21.subtract(position1);
      const Angle& angle_type =
        config.unique_type(part_type).angle(angle.type());
      ASSERT(angle_.deserialize_map().count(angle_type.model()) == 1,
        "angle model " << angle_type.model() << " not recognized.");
      en += angle_.deserialize_map()[angle_type.model()]->energy(
        relative01, relative21, angle_type);
      if (verbose_) {
        if (std::abs(en) > NEAR_ZERO) {
          const double ang = std::acos(relative01.cosine(relative21));
          const double theta0 = degrees_to_radians(angle_type.property("theta0"));
          INFO("angle " << part_index << " ijk " << angle.site(0) << " "
            << angle.site(1) << " " << angle.site(2) << " "
            << "rel01 " << relative01.str() << " "
            << "rel21 " << relative21.str() << " "
            << "ang " << ang << " "
            << "theta0 " << theta0 << " "
            << "diff " << ang - theta0);
        }
      }
    }
  }
  energy_three_body_ = en;
}

void BondVisitor::compute_four(
    const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const int part_type = part.type();
    for (const Dihedral& dihedral : config.particle_type(part_type).dihedrals()) {
      const Position& ri = part.site(dihedral.site(0)).position();
      const Position& rj = part.site(dihedral.site(1)).position();
      const Position& rk = part.site(dihedral.site(2)).position();
      const Position& rl = part.site(dihedral.site(3)).position();
      const Dihedral& dihedral_type =
        config.unique_type(part_type).dihedral(dihedral.type());
      ASSERT(dihedral_.deserialize_map().count(dihedral_type.model()) == 1,
        "dihedral model " << dihedral_type.model() << " not recognized.");
      en += dihedral_.deserialize_map()[dihedral_type.model()]->energy(
        ri, rj, rk, rl, dihedral_type);
      if (verbose_) {
        if (std::abs(en) > NEAR_ZERO) {
          FATAL("not impl");
        }
      }
    }
  }
  energy_four_body_ = en;
}

void BondVisitor::compute_all(const Select& selection,
    const Configuration& config) {
  compute_two(selection, config);
  compute_three(selection, config);
  compute_four(selection, config);
  energy_ = energy_two_body_ + energy_three_body_ + energy_four_body_;
}

void BondVisitor::compute_all(const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_all(selection, config);
}

void BondVisitor::compute_two(
    const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_two(selection, config);
}

void BondVisitor::compute_three(
    const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_three(selection, config);
}

void BondVisitor::compute_four(
    const Configuration& config,
    const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_four(selection, config);
}

}  // namespace feasst
