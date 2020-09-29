#include <string>
#include <fstream>
#include <sstream>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

Stepper::Stepper(const argtype &args) {
  args_.init(args);
  set_steps_per_write(args_.key("steps_per_write").dflt("1").integer());
  set_steps_per_update(args_.key("steps_per_update").dflt("1").integer());
  if (args_.key("file_name").used()) set_file_name(args_.str());

  if (args_.key("append").dflt("0").boolean()) {
    set_append();
  } else {
    set_no_append();
  }

  if (args_.key("clear_file").dflt("false").boolean()) {
    ASSERT(!file_name_.empty(), "file_name is a required argument with clear.");
    std::ofstream file;
    file.open(file_name_, std::ofstream::out);
    file.close();
  }
  stop_after_phase_ = args_.key("stop_after_phase").dflt("-1").integer();
  start_after_phase_ = args_.key("start_after_phase").dflt("-1").integer();
  file_name_append_phase_ =
    args_.key("file_name_append_phase").dflt("false").boolean();
  set_multistate(args_.key("multistate").dflt("0").boolean());
  is_multistate_aggregate_ =
    args_.key("multistate_aggregate").dflt("true").boolean();

  if (args_.key("num_block").used()) accumulator_.set_block(args_.integer());
  if (args_.key("num_moments").used()) accumulator_.set_moments(args_.integer());
  configuration_ = args_.key("configuration").dflt("0").integer();
}

bool Stepper::is_time(const int steps_per, int * steps_since) {
  if (steps_per > 0) {
    ++(*steps_since);
    if (*steps_since >= steps_per) {
      *steps_since = 0;
      return true;
    }
  }
  return false;
}

std::string Stepper::file_name(const Criteria& criteria) const {
  if (file_name_.empty() || !file_name_append_phase_) {
    return file_name_;
  }
  return file_name_ + "_phase" + str(criteria.phase());
}

void Stepper::printer(const std::string output, const std::string& file_name) {
  DEBUG("filename? " << file_name);
  if (file_name.empty() && !is_multistate_aggregate()) {
    std::cout << output;
  } else {
    std::ofstream file;
    if (append_) {
      file.open(file_name, std::ofstream::out | std::ofstream::app);
    } else {
      file.open(file_name, std::ofstream::out);
    }
    file << output;
    file.close();
  }
}

void Stepper::set_state(const int state) {
  state_ = state;
  if (!file_name_.empty()) {
    std::stringstream ss;
    ss << file_name_ << "_state" << state;
    file_name_ = ss.str();
  }
}

void Stepper::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  feasst_serialize_version(497, ostr);
  feasst_serialize(steps_since_update_, ostr);
  feasst_serialize(steps_since_write_, ostr);
  feasst_serialize(steps_per_update_, ostr);
  feasst_serialize(steps_per_write_, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(append_, ostr);
  feasst_serialize(stop_after_phase_, ostr);
  feasst_serialize(start_after_phase_, ostr);
  feasst_serialize(file_name_append_phase_, ostr);
  feasst_serialize(is_multistate_, ostr);
  feasst_serialize(is_multistate_aggregate_, ostr);
  feasst_serialize(state_, ostr);
  feasst_serialize(configuration_, ostr);
  feasst_serialize_fstobj(accumulator_, ostr);
  feasst_serialize_endcap("Stepper", ostr);
}

Stepper::Stepper(std::istream& istr) {
  std::string name;
  istr >> name;
  const int version = feasst_deserialize_version(istr);
  ASSERT(497 == version, version);
  feasst_deserialize(&steps_since_update_, istr);
  feasst_deserialize(&steps_since_write_, istr);
  feasst_deserialize(&steps_per_update_, istr);
  feasst_deserialize(&steps_per_write_, istr);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&append_, istr);
  feasst_deserialize(&stop_after_phase_, istr);
  feasst_deserialize(&start_after_phase_, istr);
  feasst_deserialize(&file_name_append_phase_, istr);
  feasst_deserialize(&is_multistate_, istr);
  feasst_deserialize(&is_multistate_aggregate_, istr);
  feasst_deserialize(&state_, istr);
  feasst_deserialize(&configuration_, istr);
  feasst_deserialize_fstobj(&accumulator_, istr);
  feasst_deserialize_endcap("Stepper", istr);
}

}  // namespace feasst
