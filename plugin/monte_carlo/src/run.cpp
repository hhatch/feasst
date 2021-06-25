
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

Run::Run(argtype * args) {
  num_attempts_ = integer("num_attempts", args, -1);
  until_num_particles_ = integer("until_num_particles", args, -1);
  class_name_ = "Run";
}
Run::Run(argtype args) : Run(&args) {
  check_all_used(args);
}

class MapRun {
 public:
  MapRun() {
    auto obj = MakeRun();
    obj->deserialize_map()["Run"] = obj;
  }
};

static MapRun mapper_Run = MapRun();

Run::Run(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3854, "mismatch version: " << version);
  feasst_deserialize(&num_attempts_, istr);
  feasst_deserialize(&until_num_particles_, istr);
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(num_attempts_, ostr);
  feasst_serialize(until_num_particles_, ostr);
}

void Run::perform(MonteCarlo * mc) {
  while(num_attempts_ > 0) {
    mc->attempt(1);
    --num_attempts_;
    DEBUG("num_attempts " << num_attempts_);
  }
  while(until_num_particles_ > 0 &&
        mc->configuration().num_particles() != until_num_particles_) {
    mc->attempt(1);
    DEBUG("num_particles " << mc->configuration().num_particles());
  }
}

RemoveTrial::RemoveTrial(argtype * args) {
  index_ = integer("index", args, -1);
  class_name_ = "RemoveTrial";
}
RemoveTrial::RemoveTrial(argtype args) : RemoveTrial(&args) {
  check_all_used(args);
}

class MapRemoveTrial {
 public:
  MapRemoveTrial() {
    auto obj = MakeRemoveTrial();
    obj->deserialize_map()["RemoveTrial"] = obj;
  }
};

static MapRemoveTrial mapper_RemoveTrial = MapRemoveTrial();

RemoveTrial::RemoveTrial(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3854, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
}

void RemoveTrial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(index_, ostr);
}

void RemoveTrial::perform(MonteCarlo * mc) {
  if (index_ > 0) {
    mc->remove_trial(index_);
  }
}

}  // namespace feasst
