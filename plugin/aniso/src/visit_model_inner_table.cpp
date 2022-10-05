#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

VisitModelInnerTable::VisitModelInnerTable(argtype * args) : VisitModelInner(args) {
  class_name_ = "VisitModelInnerTable";
}
VisitModelInnerTable::VisitModelInnerTable(argtype args) : VisitModelInnerTable(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void VisitModelInnerTable::precompute(Configuration * config) {
  
}

void VisitModelInnerTable::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
//  std::cout << " Info 1 part1_index " << part1_index << std::endl;
//  std::cout << " Info 1 site1_index " << site1_index << std::endl;
//  std::cout << " Info 2 part2_index " << part2_index << std::endl;
//  std::cout << " Info 2 site2_index " << site2_index << std::endl;
  TRACE("part1_index " << part1_index);
  TRACE("part2_index " << part2_index);
  TRACE("site1_index " << site1_index);
  TRACE("site2_index " << site2_index);
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
}

class MapVisitModelInnerTable {
 public:
  MapVisitModelInnerTable() {
    VisitModelInnerTable().deserialize_map()["VisitModelInnerTable"] =
      std::make_shared<VisitModelInnerTable>();
  }
};

static MapVisitModelInnerTable mapper_ = MapVisitModelInnerTable();

VisitModelInnerTable::VisitModelInnerTable(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7945, "unrecognized version: " << version);
}

void VisitModelInnerTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(7945, ostr);
}

}  // namespace feasst
