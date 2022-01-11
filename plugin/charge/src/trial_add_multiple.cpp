#include <algorithm>  // is_sorted
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/compute_add_multiple.h"

namespace feasst {

std::vector<int> ptypes(argtype * args) {
  std::vector<int> ptypes;
  int count = 0;
  std::string start("particle_type");
  std::stringstream ss;
  ss << start << count;
  while (used(ss.str(), *args)) {
    DEBUG("ss " << ss.str());
    ptypes.push_back(integer(ss.str(), args));
    ASSERT(count < 1e8, "count: " << count << " is too high");
    ++count;
    ss.str("");
    ss << start << count;
  }
  DEBUG("ptypes " << feasst_str(ptypes));
  ASSERT(std::is_sorted(ptypes.begin(), ptypes.end()),
    "ptypes not sorted: " << feasst_str(ptypes));
  return ptypes;
}

class MapTrialAddMultiple {
 public:
  MapTrialAddMultiple() {
    auto obj = MakeTrialAddMultiple();
    obj->deserialize_map()["TrialAddMultiple"] = obj;
  }
};

static MapTrialAddMultiple mapper_ = MapTrialAddMultiple();

TrialAddMultiple::TrialAddMultiple(argtype * args) : Trial(args) {
  class_name_ = "TrialAddMultiple";
  set_description("TrialAddMultiple");
  const std::vector<int> pt = ptypes(args);
  std::vector<argtype> new_args;
  set(std::make_shared<ComputeAddMultiple>(args));
  for (int p : pt) {
    argtype nag = *args;
    nag.insert({"particle_type", str(p)});
    const std::string num_steps = feasst::str("num_steps", &nag, "1");
    nag.insert({"num_steps", num_steps});
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }
  for (argtype arg : new_args) {
    add_stage(
      std::make_shared<TrialSelectParticle>(&arg),
      std::make_shared<PerturbAdd>(&arg),
      &arg);
    check_all_used(arg);
  }
}
TrialAddMultiple::TrialAddMultiple(argtype args) : TrialAddMultiple(&args) {
  //check_all_used(args);
}

TrialAddMultiple::TrialAddMultiple(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2345, "mismatch version: " << version);
}

void TrialAddMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(2345, ostr);
}

}  // namespace feasst
