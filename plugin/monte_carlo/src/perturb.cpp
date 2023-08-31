#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb.h"

namespace feasst {

Perturb::Perturb(argtype args) : Perturb(&args) { FEASST_CHECK_ALL_USED(args); }
Perturb::Perturb(argtype * args) {
  configuration_index_ = integer("configuration_index_", args, 0);
  tunable_ = Tunable(args);
}

void Perturb::before_select() {
  revert_possible_ = false;
  finalize_possible_  = false;
}

void Perturb::set_revert_possible(const bool revert_possible,
    TrialSelect * revert_select) {
  revert_possible_ = revert_possible;
  revert_select_ = revert_select;
}

void Perturb::revert(System * system) {
  if (revert_possible_) {
    FATAL("not implemented");
  }
}

void Perturb::set_finalize_possible(const bool finalize_possible,
  /// If possible, store the selection.
  TrialSelect * finalize_select) {
  finalize_possible_ = finalize_possible;
  finalize_select_ = finalize_select;
  DEBUG("finalize possible: " << finalize_possible_);
}

void Perturb::finalize(System * system) {
  //system->finalize();
  if (finalize_possible_) {
    FATAL("not implemented");
  }
}

std::string Perturb::status_header() const {
  std::stringstream ss;
  if (tunable().is_enabled()) {
    ss << ",tunable";
  }
  return ss.str();
}

std::string Perturb::status() const {
  std::stringstream ss;
  if (tunable().is_enabled()) {
    ss << "," << tunable_.value();
  }
  return ss.str();
}

std::map<std::string, std::shared_ptr<Perturb> >& Perturb::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Perturb> >* ans =
     new std::map<std::string, std::shared_ptr<Perturb> >();
  return *ans;
}

void Perturb::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Perturb> Perturb::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Perturb> Perturb::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Perturb::serialize_perturb_(std::ostream& ostr) const {
  feasst_serialize_version(903, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize_fstobj(tunable_, ostr);
}

Perturb::Perturb(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 902 && version <= 903, "mismatch version: " << version);
  if (version >= 274) {
    feasst_deserialize(&configuration_index_, istr);
  }
  feasst_deserialize_fstobj(&tunable_, istr);
}

void Perturb::perturb(
  System * system,
  TrialSelect * select,
  Random * random,
  const bool is_position_held) {
  FATAL("not implemented");
}

}  // namespace feasst
