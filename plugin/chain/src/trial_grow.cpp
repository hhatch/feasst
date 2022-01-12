#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/file.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/perturb_to_anchor.h"
#include "monte_carlo/include/perturb_distance.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/trial_compute_translate.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/select_reptate.h"
#include "chain/include/perturb_reptate.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/trial_grow.h"
#include "chain/include/select_branch.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

class MapTrialGrow {
 public:
  MapTrialGrow() {
    auto obj = MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}}});
    obj->deserialize_map()["TrialGrow"] = obj;
  }
};

static MapTrialGrow mapper_ = MapTrialGrow();

void TrialGrow::build_(std::vector<argtype> * args) {
  class_name_ = "TrialGrow";
  const int num_args = static_cast<int>(args->size());
  ASSERT(num_args > 0, "TrialGrow requires args.size: " << num_args << " > 0");
  double weight = dble("weight", &(*args)[0], -1.);
  const std::string default_num_steps = str("default_num_steps", &(*args)[0], "1");
  const std::string default_reference_index = str("default_reference_index", &(*args)[0], "-1");
  const std::string default_new_only = str("default_new_only", &(*args)[0], "false");
  // First, determine all trial types from args[0]
  std::vector<std::string> trial_types;
  if (used("transfer", (*args)[0])) {
    str("transfer", &(*args)[0]);
    trial_types.push_back("add");
    trial_types.push_back("remove");
  }
  if (used("transfer_avb", (*args)[0])) {
    str("transfer_avb", &(*args)[0]);
    trial_types.push_back("add_avb");
    trial_types.push_back("remove_avb");
  }
  if (used("regrow_avb2", (*args)[0])) {
    str("regrow_avb2", &(*args)[0]);
    trial_types.push_back("regrow_avb2_in");
    trial_types.push_back("regrow_avb2_out");
  }
  if (used("regrow", (*args)[0])) {
    str("regrow", &(*args)[0]);
    trial_types.push_back("regrow");
  }
  if (used("regrow_avb4", (*args)[0])) {
    str("regrow_avb4", &(*args)[0]);
    trial_types.push_back("regrow_avb4");
  }
  if (used("translate", (*args)[0])) {
    str("translate", &(*args)[0]);
    trial_types.push_back("translate");
  }
  if (trial_types.size() == 0) trial_types.push_back("partial_regrow");
  const std::string site = str("site", &(*args)[0], "-1");

  // Second, determine the particle_type and first site from args[0]
  const std::string particle_type = str("particle_type", &(*args)[0]);

  // Finally, add each trial to the factory
  for (const std::string trial_type : trial_types) {
    DEBUG("trial_type: " << trial_type);
    std::shared_ptr<Trial> trial = MakeTrial();
    trial->set_description("TrialGrow" + trial_type);
    if (weight > 0) {
      trial->set_weight(weight);
      //trial->set_weight(weight/static_cast<int>(trial_types.size()));
    }
    if (trial_type == "add" || trial_type == "remove" ||
        trial_type == "add_avb" || trial_type == "remove_avb" ||
        trial_type == "regrow_avb2_in" || trial_type == "regrow_avb2_out") {
      trial->set_weight(trial->weight()/2.);
    }
    std::shared_ptr<TrialCompute> compute;
    for (int iarg = 0; iarg < num_args; ++iarg) {
      DEBUG("iarg: " << iarg);
      argtype iargs = (*args)[iarg];
      std::shared_ptr<TrialSelect> select;
      std::shared_ptr<Perturb> perturb;
      if (iarg == 0 && trial_type != "partial_regrow") {
        argtype sel_args = {{"particle_type", particle_type}, {"site", site}};
        select = MakeTrialSelectParticle(sel_args);
        if (trial_type == "translate") {
          perturb = std::make_shared<PerturbTranslate>(&iargs);
          compute = std::make_shared<TrialComputeTranslate>(&iargs);
        } else if (trial_type == "add") {
          perturb = MakePerturbAdd();
          compute = MakeTrialComputeAdd();
        } else if (trial_type == "remove") {
          perturb = MakePerturbRemove();
          compute = MakeTrialComputeRemove();
        } else if (trial_type == "regrow") {
          perturb = MakePerturbAnywhere();
          compute = MakeTrialComputeMove();
        } else if (trial_type == "add_avb") {
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          auto add_avb = std::make_shared<PerturbAddAVB>(&iargs);
          ASSERT(add_avb->delay_add(), "ComputeAddAVB assumes delay_add");
          perturb = add_avb;
          compute = std::make_shared<ComputeAddAVB>();
        } else if (trial_type == "remove_avb") {
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbRemove>(
            MakePerturbMoveAVB({{{"inside", "true"}}}));
          compute = std::make_shared<ComputeRemoveAVB>();
        } else if (trial_type == "regrow_avb2_in" ||
                   trial_type == "regrow_avb2_out") {
          //argtype args_sel, args_mv, args_comp;
          //argtype args_avb2(args[iarg]);
          bool out_to_in;
          if (trial_type == "regrow_avb2_in") {
            out_to_in = true;
            iargs.insert({"out_to_in", "true"});
          } else {
            out_to_in = false;
            iargs.insert({"out_to_in", "false"});
          }
          argtype perturb_args;
          gen_avb2_args_(out_to_in, &iargs, &perturb_args);
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&perturb_args);
          compute = std::make_shared<ComputeAVB2>(&iargs);
        } else if (trial_type == "regrow_avb4") {
          //argtype args_sel, args_mv;
          gen_avb4_args_(&iargs);
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&iargs);
          compute = std::make_shared<ComputeAVB4>();
        } else {
          FATAL("unreocgnized trial_type: " << trial_type);
        }
      } else {
        int used = 0;
        if (boolean("bond", &iargs, false)) {
          ++used;
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)}});
          perturb = std::make_shared<PerturbDistance>(&iargs);
        }
        if (boolean("angle", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectAngle({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = std::make_shared<PerturbDistanceAngle>(&iargs);
        }
        if (boolean("dihedral", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectDihedral({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)},
            {"anchor_site3", str("anchor_site3", &iargs)}});
          perturb = std::make_shared<PerturbDihedral>(&iargs);
        }
        if (boolean("branch", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeSelectBranch({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"mobile_site2", str("mobile_site2", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = std::make_shared<PerturbBranch>(&iargs);
        }
        if (boolean("reptate", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)}});
          perturb = std::make_shared<PerturbToAnchor>(&iargs);
        }
        ASSERT(used == 1, "args: " << str(iargs) <<
          ". Requires one of bond, angle, dihedral, branch, reptate, etc");
        if (!compute) {
          compute = std::make_shared<TrialComputeMove>(&iargs);
        }
      }
      const std::string num_steps = str("num_steps", &iargs, default_num_steps);
      //ASSERT(trial_type != "translate" || num_steps == "1",
      //  "For " << trial_type << ", num_steps must be 1");
      argtype stage_args = {{"num_steps", num_steps},
        {"reference_index", str("reference_index", &iargs, default_reference_index)},
        {"new_only", str("new_only", &iargs, default_new_only)},
      };
      check_all_used(iargs);
      trial->add_stage(select, perturb, &stage_args);
      check_all_used(stage_args);
    }
    trial->set(compute);
    DEBUG("compute " << compute->class_name());
    add(trial);
  }
}

TrialGrow::TrialGrow(std::vector<argtype> args) : TrialFactoryNamed() {
  build_(&args);
}

class MapTrialGrowFile {
 public:
  MapTrialGrowFile() {
    auto obj = std::make_shared<TrialGrowFile>();
    obj->deserialize_map()["TrialGrowFile"] = obj;
  }
};

static MapTrialGrowFile mapper_trial_grow_file_ = MapTrialGrowFile();

void TrialGrowFile::add_(const argtype add_args, std::vector<argtype> * args) {
  argtype * arg0 = &(*args)[0];
  arg0->insert(add_args.begin(), add_args.end());
  build_(args);
  args->clear();
}

// Convert file contents to a list of argtypes to use in the above constructor
TrialGrowFile::TrialGrowFile(argtype * args) : TrialGrow() {
  const std::string file_name = str("file_name", args);
  std::vector<argtype> reformated;
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find " << file_name);
  // until end of file, find TrialGrowFile
  // check for optional empty line
  // read each line as all arguments in one stage
  // when empty line or end of file is reached, add Trial and repeat
  const bool is_found = find("TrialGrowFile", file);
  ASSERT(is_found, "TrialGrowFile not found in " << file_name);
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) {
      if (reformated.size() > 0) {
        add_(*args, &reformated);
      }
    } else {
      if (line[0] != '#') {
        reformated.push_back(line_to_argtype(line));
      }
    }
  }
  if (reformated.size() > 0) {
    add_(*args, &reformated);
  }
}
TrialGrowFile::TrialGrowFile(argtype args) : TrialGrowFile(&args) {
}
}  // namespace feasst
