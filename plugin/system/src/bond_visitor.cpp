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
  //INFO(selection.str());
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const Particle& unique_part = config.unique_type(part.type());
    const Particle& part_type = config.particle_type(part.type());
    for (int site0_index : selection.site_indices(select_index)) {
      DEBUG("site0_index " << site0_index);
      const Site& site0 = part.site(site0_index);
      for (int site1_index : part_type.bond_neighbors(site0_index)) {
        DEBUG("site1_index " << site1_index);
        const Site& site1 = part.site(site1_index);
        if (site1.is_physical()) {
          if (site0_index < site1_index ||
              !find_in_list(site1_index,
                            selection.site_indices(select_index))) {
            const Position& position0 = site0.position();
            const Position& position1 = site1.position();
            Position relative = position0;
            relative.subtract(position1);
            const Bond& bond_type = part_type.bond(site0_index, site1_index);
            const Bond& bond = unique_part.bond(bond_type.type());
            ASSERT(bond_.deserialize_map().count(bond.model()) == 1,
              "bond model " << bond.model() << " not recognized.");
            en += bond_.deserialize_map()[bond.model()]->energy(
              relative, bond);
            if (verbose_) {
              if (std::abs(en) > NEAR_ZERO) {
                INFO("bond ij " << part_index << " " << bond.site(0) << " "
                  << bond.site(1) << " sq " << relative.squared_distance());
              }
            }
          }
        }
      }
    }
  }
  energy_two_body_ = en;
}

void BondVisitor::compute_three(
    const Select& selection,
    const Configuration& config) {
  DEBUG("BondThreeBody of " << selection.str());
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const Particle& unique_part = config.unique_type(part.type());
    const Particle& part_type = config.particle_type(part.type());
    for (int site0_index : selection.site_indices(select_index)) {
      const Site& site0 = part.site(site0_index);
      for (const std::vector<int>& ang : part_type.angle_neighbors(site0_index)) {
        const int site1_index = ang[0];
        const int site2_index = ang[1];
        const Site& site1 = part.site(site1_index);
        const Site& site2 = part.site(site2_index);
        if (site1.is_physical() && site2.is_physical()) {
          if (site0_index < site1_index ||
              !find_in_list(site1_index, selection.site_indices(select_index))) {
//          if (true) {
//            if ( (site0_index < site2_index && site1_index < site2_index) ||
//                 !find_in_list(site2_index, selection.site_indices(select_index))) {
            if (true) {
              DEBUG("sites " << site0_index << " " << site1_index << " " << site2_index);
              const Position& position0 = site0.position();
              const Position& position1 = site1.position();
              const Position& position2 = site2.position();
              DEBUG("pos0 " << position0.str());
              DEBUG("pos1 " << position1.str());
              DEBUG("pos2 " << position2.str());
              Position relative01 = position0;
              Position relative21 = position2;
              relative01.subtract(position1);
              relative21.subtract(position1);
              const Angle& angle_type = part_type.angle(site0_index,
                site1_index, site2_index);
              const Angle& angle = unique_part.angle(angle_type.type());
              ASSERT(angle_.deserialize_map().count(angle.model()) == 1,
                "angle model " << angle.model() << " not recognized.");
              DEBUG(relative01.str());
              DEBUG(relative21.str());
              en += angle_.deserialize_map()[angle.model()]->energy(
                relative01, relative21, angle);
              if (verbose_) {
                if (std::abs(en) > NEAR_ZERO) {
                  const double ang = std::acos(relative01.cosine(relative21));
                  const double theta0 = degrees_to_radians(angle.property("theta0"));
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
    const Particle& unique_part = config.unique_type(part.type());
    const Particle& part_type = config.particle_type(part.type());
    for (int site0_index : selection.site_indices(select_index)) {
      const Site& site0 = part.site(site0_index);
      for (const std::vector<int>& dih : part_type.dihedral_neighbors(site0_index)) {
        const int site1_index = dih[0];
        const int site2_index = dih[1];
        const int site3_index = dih[2];
        const Site& site1 = part.site(site1_index);
        const Site& site2 = part.site(site2_index);
        const Site& site3 = part.site(site3_index);
        if (site1.is_physical() && site2.is_physical() && site3.is_physical()) {
          if (site0_index < site1_index ||
              !find_in_list(site1_index, selection.site_indices(select_index))) {
            const Position& ri = site0.position();
            const Position& rj = site1.position();
            const Position& rk = site2.position();
            const Position& rl = site3.position();
            const Dihedral& dihedral_type = part_type.dihedral(site0_index, site1_index,
              site2_index, site3_index);
            TRACE("type of dihedral " << dihedral_type.type());
            const Dihedral& dihedral = unique_part.dihedral(dihedral_type.type());
            TRACE("model of dihedral " << dihedral.model());
            ASSERT(dihedral_.deserialize_map().count(dihedral.model()) == 1,
              "dihedral model " << dihedral.model() << " not recognized.");
            en += dihedral_.deserialize_map()[dihedral.model()]->energy(
              ri, rj, rk, rl, dihedral);
            if (verbose_) {
              if (std::abs(en) > NEAR_ZERO) {
                FATAL("not impl");
              }
            }
          }
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
  //INFO("group_index " << group_index);
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
