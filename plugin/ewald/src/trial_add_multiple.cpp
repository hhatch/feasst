#include <algorithm>  // is_sorted
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/compute_add_multiple.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddMultiple(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialAddMultiple");
  // Create one stage per particle
  const std::vector<int> pt = ptypes(&args);
  std::vector<argtype> new_args;
  for (int p : pt) {
    argtype nag = args;
    nag.insert({"particle_type", str(p)});
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }
  for (argtype arg : new_args) {
    trial->add_stage(
      std::make_shared<TrialSelectParticle>(&arg),
      std::make_shared<PerturbAdd>(&arg),
      &arg);
    check_all_used(arg);
  }
  trial->set(std::make_shared<ComputeAddMultiple>());
  return trial;
}

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

}  // namespace feasst
