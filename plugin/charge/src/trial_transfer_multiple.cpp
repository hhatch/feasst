#include "utils/include/serialize.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/trial_remove_multiple.h"
#include "charge/include/trial_transfer_multiple.h"

namespace feasst {

class MapTrialTransferMultiple {
 public:
  MapTrialTransferMultiple() {
    auto obj = MakeTrialTransferMultiple();
    obj->deserialize_map()["TrialTransferMultiple"] = obj;
  }
};

static MapTrialTransferMultiple mapper_ = MapTrialTransferMultiple();

TrialTransferMultiple::TrialTransferMultiple(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransferMultiple";
  argtype rm_args = *args, add_args = *args;
  str("shift", &add_args, ""); // remove shift
  add(MakeTrialAddMultiple(add_args));
  add(MakeTrialRemoveMultiple(rm_args));
}
TrialTransferMultiple::TrialTransferMultiple(argtype args) : TrialTransferMultiple(&args) {
  //check_all_used(args);
}

}  // namespace feasst
