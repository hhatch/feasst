#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/potential.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "confinement/include/zero_potential.h"

namespace feasst {

ZeroPotential::ZeroPotential(argtype * args) {
  class_name_ = "ZeroPotential";
  configuration_index_ = integer("configuration_index", args, 0);
}
ZeroPotential::ZeroPotential(argtype args) : ZeroPotential(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(ZeroPotential,);

ZeroPotential::ZeroPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1497, "mismatch version: " << version);
  feasst_deserialize(&configuration_index_, istr);
}

void ZeroPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(1497, ostr);
  feasst_serialize(configuration_index_, ostr);
}

void ZeroPotential::run(MonteCarlo * mc) {
  const double current_energy = mc->criteria().current_energy();
  mc->add(MakePotential(
  _to_reference(MakePotential(args_),
                       reference_index_,
                       configuration_index_);
}

}  // namespace feasst
