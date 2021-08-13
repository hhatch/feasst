
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
  class_name_ = "RemoveTrial";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
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
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveTrial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveTrial::perform(MonteCarlo * mc) {
  if (!name_.empty()) {
    for (int trial = 0; trial < mc->trials().num(); ++trial) {
      if (mc->trial(trial).class_name() == name_) {
        ASSERT(index_ < 0 || trial == index_,
          "RemoveTrial cannot specify both index and name");
        index_ = trial;
        break;
      }
    }
  }
  if (index_ >= 0) {
    mc->remove_trial(index_);
  }
  if (all_) {
    for (int i = 0; mc->trials().num() > 0; ++i) {
      mc->remove_trial(0);
      ASSERT(i < 1e5, "too many trials. Infinite loop?");
    }
  }
}

RemoveModify::RemoveModify(argtype * args) {
  class_name_ = "RemoveModify";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveModify::RemoveModify(argtype args) : RemoveModify(&args) {
  check_all_used(args);
}

class MapRemoveModify {
 public:
  MapRemoveModify() {
    auto obj = MakeRemoveModify();
    obj->deserialize_map()["RemoveModify"] = obj;
  }
};

static MapRemoveModify mapper_RemoveModify = MapRemoveModify();

RemoveModify::RemoveModify(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2045, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveModify::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2045, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveModify::perform(MonteCarlo * mc) {
  DEBUG("name " << name_);
  if (!name_.empty()) {
    for (int modify = 0; modify < mc->num_modifiers(); ++modify) {
      DEBUG("mod " << mc->modify(modify).class_name());
      if (mc->modify(modify).class_name() == name_) {
        ASSERT(index_ < 0 || modify == index_,
          "RemoveModify cannot specify both index and name");
        index_ = modify;
        DEBUG("removing " << modify);
        break;
      }
    }
  }
  DEBUG("index " << index_);
  if (index_ >= 0) {
    DEBUG("removing " << index_);
    mc->remove_modify(index_);
  }
  if (all_) {
    for (int i = 0; mc->num_modifiers() > 0; ++i) {
      mc->remove_modify(0);
      ASSERT(i < 1e5, "Infinite loop?");
    }
  }
}

WriteCheckpoint::WriteCheckpoint(argtype * args) {
  class_name_ = "WriteCheckpoint";
}
WriteCheckpoint::WriteCheckpoint(argtype args) : WriteCheckpoint(&args) {
  check_all_used(args);
}

class MapWriteCheckpoint {
 public:
  MapWriteCheckpoint() {
    auto obj = MakeWriteCheckpoint();
    obj->deserialize_map()["WriteCheckpoint"] = obj;
  }
};

static MapWriteCheckpoint mapper_WriteCheckpoint = MapWriteCheckpoint();

WriteCheckpoint::WriteCheckpoint(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5694, "mismatch version: " << version);
}

void WriteCheckpoint::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(5694, ostr);
}

void WriteCheckpoint::perform(MonteCarlo * mc) {
  mc->write_checkpoint();
}

}  // namespace feasst
