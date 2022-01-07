#include <sstream>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_cell.h"
#include "system/include/model_two_body.h"
#include "system/include/model_one_body.h"

namespace feasst {

VisitModelCell::VisitModelCell(argtype * args) : VisitModel() {
  class_name_ = "VisitModelCell";
  DEBUG("parse cells");
  if (used("min_length", *args)) {
    min_length_ = dble("min_length", args);
    group_index_ = integer("cell_group", args, 0);
    ASSERT(group_index_ >= 0, "invalid group_index: " << group_index_);
  }
}
VisitModelCell::VisitModelCell(argtype args) : VisitModelCell(&args) {
  check_all_used(args);
}
VisitModelCell::VisitModelCell(std::shared_ptr<VisitModelInner> inner,
  argtype args) : VisitModelCell(args) {
  set_inner(inner);
}

int VisitModelCell::cell_id(const Domain& domain,
                            const Position& position) const {
  Position pos = position;
  domain.wrap(&pos);
  pos.divide(domain.side_lengths());
  return cells_.id(pos.coord());
}

// HWH note if there are problems with scaled coordinates here, it probably
// means there is an issue with wrapping. As currently implemented, translations
// automatically wrap. So if you're doing a test without them you might run
// into this issue.
int VisitModelCell::cell_id_opt_(const Domain& domain,
                                 const Position& position) {
//  Position scaled(position);
//  DEBUG("scaled before wrap " << scaled.str() << " pos " << position.str() <<
//    " box " << side_lengths().str());
  domain.wrap_opt(position, opt_origin_, &opt_rel_, &opt_pbc_, &opt_r2_);
  //wrap(&scaled);
  DEBUG("opt_rel_ after wrap " << opt_rel_.str() << " pos " << position.str());
  opt_rel_.divide(domain.side_lengths());
  DEBUG("opt_rel_ " << opt_rel_.str() << " pos " << position.str());
  return cells_.id(opt_rel_.coord());
}

void VisitModelCell::precompute(Configuration * config) {
  VisitModel::precompute(config);
  ASSERT(config->domain().side_lengths().size() > 0,
    "cannot define cells before domain sides");
  ASSERT(!config->domain().is_tilted(), "implement triclinic");
  if (cells_.type() == -1) {
    Cells cells;
    cells.create(min_length_, config->domain().side_lengths().coord());
    cells.set_type(config->num_cell_lists());
    config->increment_num_cell_lists();
    cells.set_group(group_index_);
    if (cells.num_total() > 0) {
      cells_ = cells;
    } else {
      FATAL("Requested cell list rejected: min_length:" << min_length_ <<
            " did not meet requirements.");
    }
    opt_origin_.set_to_origin(config->dimension());
    opt_rel_.set_to_origin(config->dimension());
    opt_pbc_.set_to_origin(config->dimension());
    position_tracker_(config->group_selects()[group_index_], config);
  }
  check(*config);
}

void VisitModelCell::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(group_index == group_index_, "not equivalent");
  init_relative_(domain, &relative_, &pbc_);

  /*
    Loop index nomenclature
    ends in 1 or 2 to represent the pair
    cell -> id of cell
    select -> selection inside each cell
    select_index -> index of selection
    part -> particle
    part_index -> index of particle
    site -> site
   */

  // loop through neighboring cells where cell1 < cell2 only
  for (int cell1 = 0; cell1 < cells_.num_total(); ++cell1) {
    const Select& select1 = cells_.particles()[cell1];
    for (int cell2 : cells_.neighbor()[cell1]) {
      if (cell1 < cell2) {
        const Select& select2 = cells_.particles()[cell2];
        for (int select1_index = 0;
             select1_index < select1.num_particles();
             ++select1_index) {
          const int part1_index = select1.particle_index(select1_index);
          for (int select2_index = 0;
               select2_index < select2.num_particles();
               ++select2_index) {
            const int part2_index = select2.particle_index(select2_index);
            if (part1_index != part2_index) {
              for (int site1_index : select1.site_indices(select1_index)) {
                for (int site2_index : select2.site_indices(select2_index)) {
                  get_inner_()->compute(part1_index, site1_index, part2_index,
                                        site2_index, config, model_params,
                                        model, false, &relative_, &pbc_);
                }
              }
            }
          }
        }
      }
    }
  }

  // loop through the same cell only
  for (int cell1 = 0; cell1 < cells_.num_total(); ++cell1) {
    const Select& select = cells_.particles()[cell1];
    for (int select1_index = 0;
         select1_index < select.num_particles() - 1;
         ++select1_index) {
      const int part1_index = select.particle_index(select1_index);
      for (int select2_index = select1_index + 1;
           select2_index < select.num_particles();
           ++select2_index) {
        const int part2_index = select.particle_index(select2_index);
        if (part1_index != part2_index) {
          for (int site1_index : select.site_indices(select1_index)) {
            for (int site2_index : select.site_indices(select2_index)) {
              get_inner_()->compute(part1_index, site1_index, part2_index,
                                    site2_index, config, model_params, model,
                                    false, &relative_, &pbc_);
            }
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

void VisitModelCell::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("visiting model");
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(selection.num_particles() == 1, "assumes 1 particle selection");
  ASSERT(group_index == group_index_, "not equivalent");
  init_relative_(domain, &relative_, &pbc_);
  for (int select1_index = 0;
       select1_index < selection.num_particles();
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    const Particle& part1 = config->select_particle(part1_index);
    for (int site1_index : selection.site_indices(select1_index)) {
      const Site& site1 = part1.site(site1_index);
      const int cell1_index = cell_id_opt_(domain, site1.position());
      for (int cell2_index : cells_.neighbor()[cell1_index]) {
        const Select& cell2_parts = cells_.particles()[cell2_index];
        for (int select2_index = 0;
             select2_index < cell2_parts.num_particles();
             ++select2_index) {
          const int part2_index = cell2_parts.particle_index(select2_index);
          if (part1_index != part2_index) {
            TRACE("indices " <<
                  feasst_str(cell2_parts.site_indices(select2_index)));
            for (int site2_index : cell2_parts.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                   site1_index << " " << site2_index);
              get_inner_()->compute(part1_index, site1_index, part2_index,
                                    site2_index, config, model_params, model,
                                    false, &relative_, &pbc_);
            }
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

void VisitModelCell::position_tracker_(const Select& select,
    Configuration * config) {
  for (int spindex = 0; spindex < select.num_particles(); ++spindex) {
    const int particle_index = select.particle_index(spindex);
    for (const int site_index : select.site_indices(spindex)) {
      ASSERT(site_index >= 0, "index error");
      DEBUG("update cells");
      DEBUG("group " << cells_.group());
      const int group_index = cells_.group();
      ASSERT(group_index >= 0, "error");
      DEBUG("config " << config);
      DEBUG("particle_index " << particle_index);
      DEBUG("num particles " << config->num_particles());
      const Particle& part = config->select_particle(particle_index);
      ASSERT(config, "error");
      DEBUG("group_index: " << group_index);
      DEBUG("group selects size: " << config->group_selects().size());
      const Select& sel = config->group_selects()[group_index];
      DEBUG("sel " << sel.str());
      DEBUG("is group empty: " << sel.is_group_empty());
      const Group& group = config->group_selects()[group_index].group();
      if (group.is_in(part)) {
        const Site& site = part.site(site_index);
        if (group.is_in(site)) {
          const int cell_new = cell_id_opt_(config->domain(), site.position());
          if (one_site_select_.num_particles() == 0) {
            one_site_select_.add_site(0, 0);
          }
          one_site_select_.set_particle(0, particle_index);
          one_site_select_.set_site(0, 0, site_index);
          ParticleFactory * particles = config->get_particles_();
          Site * sitep = particles->get_particle(particle_index)->get_site(site_index);
          if (cells_.type() < site.num_cells()) {
            DEBUG(cells_.type());
            const int cell_old = site.cell(cells_.type());
            DEBUG("index " << particle_index << " " << site_index);
            DEBUG("new cell " << cell_new << " old cell " << cell_old);
//            DEBUG("before new cell set: " <<
//              particles->particle(particle_index).site(
//              site_index).property("cell0"));
            sitep->set_cell(cells_.type(), cell_new);
            DEBUG(one_site_select_.str());
            DEBUG(cells_.num_total());
            cells_.update(one_site_select_, cell_new, cell_old);
          } else {
            sitep->add_cell(cell_new);
            DEBUG("adding to cell list cllnw "
              << cell_new << " si " << site_index);
            cells_.add(one_site_select_, cell_new);
          }
        }
      }
    }
  }
}

void VisitModelCell::finalize(const Select& select, Configuration * config) {
  VisitModel::finalize(select, config);
  if (select.trial_state() == 2) {
    // remove particles from cell
    for (const int particle_index : select.particle_indices()) {
      // note: somewhat derivative of position_tracker
      const int group_index = cells_.group();
      const Particle& part = config->select_particle(particle_index);
      const Group& group = config->group_selects()[group_index].group();
      if (group.is_in(part)) {
        for (int site_index = 0; site_index < part.num_sites(); ++site_index) {
          const Site& site = part.site(site_index);
          if (group.is_in(site)) {
            if (cells_.type() < site.num_cells()) {
              const int cell_old = site.cell(cells_.type());
              Select select;
              select.add_site(particle_index, site_index);
              cells_.remove(select, cell_old);
            }
          }
        }
      }
    }
  } else {
    position_tracker_(select, config);
  }
}

void VisitModelCell::check(const Configuration& config) const {
  VisitModel::check(config);
  // for each site in config that has cells, check that the cell is correct.
  int num_sites_in_cell = 0;
  for (const int part_index : config.selection_of_all().particle_indices()) {
    for (const Site& site : config.select_particle(part_index).sites()) {
      if (site.num_cells() > cells_.type()) {
        ++num_sites_in_cell;
        const int old_cell = site.cell(cells_.type());
        const int cur_cell = cell_id(config.domain(), site.position());
        ASSERT(old_cell == cur_cell,
          "old_cell: " << old_cell << " != cur_cell: " << cur_cell);
      }
    }
  }
  ASSERT(num_sites_in_cell == cells_.num_sites(),
    "num sites with cells: " << num_sites_in_cell << " != " <<
    cells_.num_sites());
}

class MapVisitModelCell {
 public:
  MapVisitModelCell() {
    auto obj = MakeVisitModelCell({{"min_length", "1"}});
    obj->deserialize_map()["VisitModelCell"] = obj;
  }
};

static MapVisitModelCell mapper_ = MapVisitModelCell();

VisitModelCell::VisitModelCell(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(755 == version, version);
  feasst_deserialize(&min_length_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize_fstobj(&opt_origin_, istr);
  feasst_deserialize_fstobj(&opt_rel_, istr);
  feasst_deserialize_fstobj(&opt_pbc_, istr);
  feasst_deserialize_fstobj(&cells_, istr);
}

void VisitModelCell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(755, ostr);
  feasst_serialize(min_length_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstobj(opt_origin_, ostr);
  feasst_serialize_fstobj(opt_rel_, ostr);
  feasst_serialize_fstobj(opt_pbc_, ostr);
  feasst_serialize_fstobj(cells_, ostr);
}

}  // namespace feasst
