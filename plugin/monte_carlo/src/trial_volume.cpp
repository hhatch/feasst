#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_volume.h"

namespace feasst {

class MapTrialVolume {
 public:
  MapTrialVolume() {
    auto obj = MakeTrialVolume();
    obj->deserialize_map()["TrialVolume"] = obj;
  }
};

static MapTrialVolume mapper_ = MapTrialVolume();

TrialVolume::TrialVolume(argtype * args) : Trial(args) {
  class_name_ = "TrialVolume";
  set_description("TrialVolume");
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbVolume>(args),
    args);
  set(MakeTrialComputeVolume());
}
TrialVolume::TrialVolume(argtype args) : TrialVolume(&args) {
  check_all_used(args);
}

TrialVolume::TrialVolume(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8945, "mismatch version: " << version);
}

void TrialVolume::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(8945, ostr);
}

}  // namespace feasst
