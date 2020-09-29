#include "steppers/include/extensive_moments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"

namespace feasst {

class MapExtensiveMoments {
 public:
  MapExtensiveMoments() {
    auto obj = MakeExtensiveMoments();
    obj->deserialize_map()["ExtensiveMoments"] = obj;
  }
};

static MapExtensiveMoments mapper_ = MapExtensiveMoments();

ExtensiveMoments::ExtensiveMoments(const argtype &args) : Analyze(args) {
  args_.init(args);
  max_order_ = args_.key("max_order").dflt("3").integer();
}

void ExtensiveMoments::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const int num_ptypes = system->configuration().num_particle_types();
  resize(max_order_ + 1,
         max_order_ + 1,
         num_ptypes,
         max_order_ + 1,
         num_ptypes,
         &moments_);
  u_p_.resize(max_order_ + 1);
  resize(num_ptypes, max_order_ + 1, &n_i_);
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string ExtensiveMoments::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  //ss << "max_order," << max_order_ << std::endl;
  return ss.str();
}

void ExtensiveMoments::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const int num_ptypes = system.configuration().num_particle_types();

  // precompute powers
  const double energy = criteria.current_energy();
  u_p_[0] = 1.;
  n_i_[0] = std::vector<double>(max_order_ + 1, 1.);
  for (int order = 1; order <= max_order_; ++order) {
    u_p_[order] = energy*u_p_[order - 1];
    for (int ptype = 0; ptype < num_ptypes; ++ptype) {
      n_i_[ptype][order] = system.configuration().num_particles_of_type(ptype)*
                           n_i_[ptype][order - 1];
    }
  }

  // accumulate moments
  for (int p = 0; p <= max_order_; ++p) {
  for (int m = 0; m <= max_order_; ++m) {
  for (int k = 0; k < num_ptypes; ++k) {
  for (int j = 0; j <= max_order_; ++j) {
  for (int i = 0; i < num_ptypes; ++i) {
    moments_[p][m][k][j][i].accumulate(n_i_[i][j]*n_i_[k][m]*u_p_[p]);
  }}}}}
}

std::string ExtensiveMoments::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  feasst_serialize_fstobj(moments_, ss);
  return ss.str();
}

void ExtensiveMoments::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1647, ostr);
  feasst_serialize(max_order_, ostr);
  feasst_serialize_fstobj(moments_, ostr);
  feasst_serialize(u_p_, ostr);
  feasst_serialize(n_i_, ostr);
}

ExtensiveMoments::ExtensiveMoments(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1647, "mismatch version:" << version);
  feasst_deserialize(&max_order_, istr);
  feasst_deserialize_fstobj(&moments_, istr);
  feasst_deserialize(&u_p_, istr);
  feasst_deserialize(&n_i_, istr);
}

ExtensiveMoments::ExtensiveMoments(const Analyze& extensive_moments) {
  std::stringstream ss;
  extensive_moments.serialize(ss);
  *this = ExtensiveMoments(ss);
}

}  // namespace feasst
