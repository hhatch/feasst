#include <cmath>
#include <vector>
#include "utils/include/serialize.h"
#include "utils/include/utils.h"  // find_in_list
#include "configuration/include/configuration.h"
#include "system/include/visit_model_intra.h"
#include "system/include/model_two_body.h"

namespace feasst {

VisitModelIntra::VisitModelIntra(argtype * args) : VisitModel() {
  class_name_ = "VisitModelIntra";
  set_cutoff(integer("cutoff", args, -1));
}
VisitModelIntra::VisitModelIntra(argtype args) : VisitModelIntra(&args) {
  check_all_used(args);
}

void VisitModelIntra::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  TRACE("intra particle energy_of_selection");
  ASSERT(group_index == 0,
    "need to implement site1 loop filtering particles by group");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  for (int sp1index = 0;
       sp1index < static_cast<int>(selection.particle_indices().size());
       ++sp1index) {
    const int part1_index = selection.particle_index(sp1index);
    TRACE("particle: " << part1_index);
    const Particle& part1 = config->select_particle(part1_index);
    // the first site loop is over all sites in part1 and group_index
    // the second is all sites in selection
    // HWH optimize this

    // here we use excluded to account for chain regrowth, etc.
    // exclude the particles which haven't been grown yet.
    // or exclude particles which will form new bonds (reptate).
    Select sites1;
    sites1.add_sites(selection.particle_index(sp1index),
                     selection.site_indices(sp1index));
    if (selection.excluded()) {
      sites1.remove(*(selection.excluded()));
      TRACE("excluded " << selection.excluded()->str());
    }
    TRACE("sites1: " << sites1.str());
    const std::vector<int>& site1_indices = sites1.site_indices(0);

    Select sites2;
    sites2.add_particle(part1, part1_index);
    if (selection.excluded()) {
      sites2.remove(*(selection.excluded()));
      TRACE("excluded " << selection.excluded()->str());
    }
    TRACE("sites2: " << sites2.str());
    const std::vector<int>& site2_indices = sites2.site_indices(0);
    for (const int site1_index : site1_indices) {
      for (const int site2_index : site2_indices) {
        // if sites in particle selection > 1, attempt the following check.
        // if site2 is in selection, then require site1 < site2
        if (!find_in_list(site2_index, site1_indices) ||
            site1_index < site2_index) {
          bool include = false;
          if (selection.old_bond()) {
            if (site2_index ==
                selection.old_bond()->site_indices()[sp1index][0]) {
              include = true;
            }
          }
          bool exclude = false;
          if (selection.new_bond()) {
            if (site2_index ==
                selection.new_bond()->site_indices()[sp1index][0]) {
              exclude = true;
            }
          }

          // forced exclude takes precedent over forced include
          if ( (include || std::abs(site1_index - site2_index) > cutoff_) &&
               (!exclude) ) {
            TRACE("sites: " << site1_index << " " << site2_index);
            get_inner_()->compute(part1_index, site1_index, part1_index,
              site2_index, config, model_params, model, false, &relative_,
              &pbc_);
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

class MapVisitModelIntra {
 public:
  MapVisitModelIntra() {
    VisitModelIntra().deserialize_map()["VisitModelIntra"] =
      std::make_shared<VisitModelIntra>();
  }
};

static MapVisitModelIntra mapper_ = MapVisitModelIntra();

VisitModelIntra::VisitModelIntra(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(754 == version, version);
  feasst_deserialize(&cutoff_, istr);
}

void VisitModelIntra::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(754, ostr);
  feasst_serialize(cutoff_, ostr);
}

void VisitModelIntra::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  compute(model, model_params, config->selection_of_all(),
          config, group_index);
}

}  // namespace feasst
