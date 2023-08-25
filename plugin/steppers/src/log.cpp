#include "steppers/include/log.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapLog {
 public:
  MapLog() {
    Log().deserialize_map()["Log"] = MakeLog();
  }
};

static MapLog mapper_ = MapLog();

Log::Log(argtype args) : Log(&args) { FEASST_CHECK_ALL_USED(args); }
Log::Log(argtype * args) : AnalyzeWriteOnly(args) {
  if (boolean("append", args, true)) {
    set_append();
  } else {
    ERROR("append is required");
  }
  max_precision_ = boolean("max_precision", args, false);
  include_bonds_ = boolean("include_bonds", args, false);
}

void Log::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string Log::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << system.status_header();
  ss << criteria.status_header(system);
  if (include_bonds_) {
    ss << ",BondTwoBody,BondThreeBody,BondFourBody";
  }
  // print number of trials here instead of TrialFactory header because
  // multiple factories makes it redundant.
  ss << ",trial"
     << trial_factory.status_header()
     << std::endl;
  return ss.str();
}

std::string Log::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
  ss << system.status();
  ss << criteria.status(max_precision_);
  if (include_bonds_) {
    bond_visitor_.compute_all(system.configuration());
    if (max_precision_) {
      ss << "," << MAX_PRECISION << bond_visitor_.energy_two_body()
         << "," << MAX_PRECISION << bond_visitor_.energy_three_body()
         << "," << MAX_PRECISION << bond_visitor_.energy_four_body();
    } else {
      ss << "," << bond_visitor_.energy_two_body()
         << "," << bond_visitor_.energy_three_body()
         << "," << bond_visitor_.energy_four_body();
    }
  }
  ss << "," << trial_factory.num_attempts()
     << trial_factory.status()
     << std::endl;
  return ss.str();
}

void Log::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(668, ostr);
  feasst_serialize(max_precision_, ostr);
  feasst_serialize(include_bonds_, ostr);
  feasst_serialize_fstobj(bond_visitor_, ostr);
}

Log::Log(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 668, "version mismatch:" << version);
  feasst_deserialize(&max_precision_, istr);
  feasst_deserialize(&include_bonds_, istr);
  feasst_deserialize_fstobj(&bond_visitor_, istr);
}

}  // namespace feasst
