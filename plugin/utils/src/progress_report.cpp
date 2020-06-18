#include <iostream>
#include <iomanip>
#include <fstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "utils/include/progress_report.h"

namespace feasst {

void ProgressReport::reset() {
  last_percent_ = -1.;
  current_ = 0.;
}

ProgressReport::ProgressReport(const argtype &args) {
  reset();
  Arguments args_(args);
  num_ = args_.key("num").dflt("1").integer();
  percent_per_write_ = args_.key("percent").dflt("0.1").dble();
  if (args_.key("file_name").used()) file_name_ = args_.str();
  ASSERT(percent_per_write_ > 0.01, "percent_per_write: " << percent_per_write_
    << " must be > 0.01");
}

void ProgressReport::serialize(std::ostream& ostr) const {
  feasst_serialize_version(5756, ostr);
  feasst_serialize(num_, ostr);
  feasst_serialize(current_, ostr);
  feasst_serialize(percent_per_write_, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(last_percent_, ostr);
}

ProgressReport::ProgressReport(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5756, "version mismatch: " << version);
  feasst_deserialize(&num_, istr);
  feasst_deserialize(&current_, istr);
  feasst_deserialize(&percent_per_write_, istr);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&last_percent_, istr);
}

double ProgressReport::percent() const {
  return static_cast<double>(current_)/static_cast<double>(num_);
}

void ProgressReport::check() {
  if (current_ == 0) {
    starting_hours_ = cpu_hours();
  }
  ++current_;
  if (current_ == num_ || percent() >= last_percent_ + percent_per_write_) {
    write();
    last_percent_ = percent();
  }
}

void ProgressReport::write() {
  std::stringstream ss;
  if (current_ == 1) {
    ss << "# Beginning progress report.";
  } else if (current_ == num_) {
    ss << "# 100 % complete.";
  } else {
    const double elapsed_hours = cpu_hours() - starting_hours_;
    DEBUG("elapsed " << elapsed_hours);
    const double percent_per_hours = percent()/elapsed_hours;
    DEBUG("percent_per_hours " << percent_per_hours);
    const double remaining_hours = (1. - percent())/percent_per_hours;
    ss << "# " << std::setprecision(3) <<  percent()*100 << "% " << elapsed_hours
       << " hours elapsed. Estimated " << remaining_hours
       << " hours remain.";
  }
  DEBUG("filename? " << file_name_);
  if (file_name_.empty()) {
    std::cout << ss.str() << std::endl;
  } else {
    std::ofstream file;
    file.open(file_name_, std::ofstream::out | std::ofstream::app);
    file << ss.str() << std::endl;
    file.close();
  }
}

}  // namespace feasst
