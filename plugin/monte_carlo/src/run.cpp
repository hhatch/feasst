
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

Run::Run(argtype * args) {
  num_attempts_ = integer("num_attempts", args, -1);
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

static MapRun mapper_ = MapRun();

Run::Run(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3854, "mismatch version: " << version);
  feasst_deserialize(&num_attempts_, istr);
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(num_attempts_, ostr);
}

void Run::perform(MonteCarlo * mc) {
  while(num_attempts_ > 0) {
    mc->attempt(1);
    --num_attempts_;
  }
}

}  // namespace feasst
