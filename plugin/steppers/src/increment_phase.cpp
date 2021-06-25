#include "utils/include/serialize.h"
#include "steppers/include/increment_phase.h"

namespace feasst {

class MapIncrementPhase {
 public:
  MapIncrementPhase() {
    IncrementPhase().deserialize_map()["IncrementPhase"] = MakeIncrementPhase();
  }
};

static MapIncrementPhase mapper_energy_check_ = MapIncrementPhase();

IncrementPhase::IncrementPhase(argtype * args) : ModifyUpdateOnly(args) {
  num_trials_ = integer("num_trials", args, -1);
}
IncrementPhase::IncrementPhase(argtype args) : IncrementPhase(&args) {
  check_all_used(args);
}

void IncrementPhase::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  if (num_trials_ != -1) {
    if (trial_factory->num_attempts() > num_trials_) {
      criteria->increment_phase();
      num_trials_ = -1;  // disable further increments
    }
  }
}

void IncrementPhase::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(4958, ostr);
  feasst_serialize(num_trials_, ostr);
}

IncrementPhase::IncrementPhase(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4958, "version mismatch: " << version);
  feasst_deserialize(&num_trials_, istr);
}

}  // namespace feasst
