#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_transfer.h"

namespace feasst {

class MapTrialTransfer {
 public:
  MapTrialTransfer() {
    auto obj = MakeTrialTransfer();
    obj->deserialize_map()["TrialTransfer"] = obj;
  }
};

static MapTrialTransfer mapper_ = MapTrialTransfer();

TrialTransfer::TrialTransfer(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransfer";
  argtype orig_args = *args;
  add(std::make_shared<TrialAdd>(args));
  add(MakeTrialRemove(orig_args));
}
TrialTransfer::TrialTransfer(argtype args) : TrialTransfer(&args) {
  check_all_used(args);
}

}  // namespace feasst