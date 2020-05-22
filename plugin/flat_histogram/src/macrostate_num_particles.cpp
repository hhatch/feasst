#include <algorithm>
#include "utils/include/serialize.h"
#include "flat_histogram/include/macrostate_num_particles.h"

namespace feasst {

MacrostateNumParticles::MacrostateNumParticles(const Histogram& histogram,
    const argtype& args) : Macrostate(histogram, args) {
  num_ = ConstrainNumParticles(
    {{"type", args_.key("particle_type").dflt("-1").str()}});
  ASSERT(num_.type() >= -1, "particle_type: " << num_.type());
}

double MacrostateNumParticles::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  return num_.num_particles(system, acceptance);
}

class MapMacrostateNumParticles {
 public:
  MapMacrostateNumParticles() {
    auto hist = MakeHistogram({{"width", "1"}, {"max", "1"}});
    MacrostateNumParticles(*hist).deserialize_map()["MacrostateNumParticles"] =
      MakeMacrostateNumParticles(*hist);
  }
};

static MapMacrostateNumParticles mapper_ = MapMacrostateNumParticles();

std::shared_ptr<Macrostate> MacrostateNumParticles::create(std::istream& istr) const {
  return std::make_shared<MacrostateNumParticles>(istr);
}

MacrostateNumParticles::MacrostateNumParticles(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 204, "version mismatch: " << version);
  feasst_deserialize_fstobj(&num_, istr);
}

void MacrostateNumParticles::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(204, ostr);
  feasst_serialize_fstobj(num_, ostr);
}

}  // namespace feasst
