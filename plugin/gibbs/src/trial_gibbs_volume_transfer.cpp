#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"
#include "gibbs/include/trial_gibbs_volume_transfer.h"

namespace feasst {

class MapTrialGibbsVolumeTransferOneWay {
 public:
  MapTrialGibbsVolumeTransferOneWay() {
    auto obj = MakeTrialGibbsVolumeTransferOneWay({{"to_configuration_index", "1"}});
    obj->deserialize_map()["TrialGibbsVolumeTransferOneWay"] = obj;
  }
};

static MapTrialGibbsVolumeTransferOneWay mapper_trial_remove_avb_ = MapTrialGibbsVolumeTransferOneWay();

TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsVolumeTransferOneWay";
  set_description("TrialGibbsVolumeTransferOneWay");
  const int to_configuration_index = integer("to_configuration_index", args);
  const int configuration_index = integer("configuration_index", args, 0);
  argtype add_args = *args;
  add_args.insert({"configuration_index", str(to_configuration_index)});
  add_stage(
    std::make_shared<TrialSelectParticle>(&add_args),
    std::make_shared<PerturbAdd>(),
    &add_args);
  FEASST_CHECK_ALL_USED(add_args);
  args->insert({"configuration_index", str(configuration_index)});
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbRemove>(),
    args);
  set(MakeComputeGibbsVolumeTransfer());
}
TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(argtype args) : TrialGibbsVolumeTransferOneWay(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9234, "mismatch version: " << version);
}

void TrialGibbsVolumeTransferOneWay::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9234, ostr);
}

class MapTrialGibbsVolumeTransfer {
 public:
  MapTrialGibbsVolumeTransfer() {
    auto obj = MakeTrialGibbsVolumeTransfer();
    obj->deserialize_map()["TrialGibbsVolumeTransfer"] = obj;
  }
};

static MapTrialGibbsVolumeTransfer mapper_trial_transfer_avb__ = MapTrialGibbsVolumeTransfer();

TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialGibbsVolumeTransfer";
  const int config0 = integer("configuration_index0", args, 0);
  const int config1 = integer("configuration_index1", args, 1);
  //INFO("args " << str(*args));
  //ASSERT(!used("configuration_index", *args),
  //  "Do not use argument:configuration_index. Use configuration_index0 or 1.");
  argtype orig_args = *args;
  orig_args.insert({"configuration_index", str(config0)});
  orig_args.insert({"to_configuration_index", str(config1)});
  args->insert({"configuration_index", str(config1)});
  args->insert({"to_configuration_index", str(config0)});
  auto trial_add = MakeTrialGibbsVolumeTransferOneWay(orig_args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = std::make_shared<TrialGibbsVolumeTransferOneWay>(args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype args) : TrialGibbsVolumeTransfer(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
