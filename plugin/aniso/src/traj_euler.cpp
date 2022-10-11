#include "utils/include/serialize.h"
#include "aniso/include/traj_euler.h"

namespace feasst {

class MapTrajEuler {
 public:
  MapTrajEuler() {
    auto obj = MakeTrajEuler({{"file_name", "place_holder"}});
    obj->deserialize_map()["TrajEuler"] = obj;
  }
};

static MapTrajEuler mapper_ = MapTrajEuler();

TrajEuler::TrajEuler(argtype * args) : AnalyzeWriteOnly(args) {
  set_append();
  ASSERT(!file_name().empty(), "file name is required");
  args->insert({"append", "true"}); // always append
  xyz_ = FileXYZEuler(args);
}
TrajEuler::TrajEuler(argtype args) : TrajEuler(&args) { FEASST_CHECK_ALL_USED(args); }

void TrajEuler::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const std::string name = file_name(*criteria);
  ASSERT(!name.empty(), "file name required. Did you forget to " <<
    "Analyze::set_file_name()?");

  // write xyz
  if (state() == criteria->state()) {
    xyz_.write(name, system->configuration());
  }
}

std::string TrajEuler::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  xyz_.write(file_name(criteria), system.configuration());
  return std::string("");
}

void TrajEuler::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(8246, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
}

TrajEuler::TrajEuler(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8246, "version mismatch:" << version);
  feasst_deserialize_fstobj(&xyz_, istr);
}

}  // namespace feasst
