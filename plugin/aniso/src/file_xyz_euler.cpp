
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "aniso/include/file_xyz_euler.h"

namespace feasst {

FileXYZEuler::FileXYZEuler(argtype * args) {
  group_index_ = integer("group_index", args, 0);
  append_ = boolean("append", args, false);
}
FileXYZEuler::FileXYZEuler(argtype args) : FileXYZEuler(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void FileXYZEuler::load(const std::string file_name, Configuration * config) const {
  FATAL("not implemented");
}

PrinterXYZEuler::PrinterXYZEuler(std::shared_ptr<std::ofstream> file,
    const int num_places) : file_(file) {
  num_places_ = num_places;
}

void PrinterXYZEuler::work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) {
  (*file_.get()) << site.type() << " ";
  (*file_.get()) << std::setprecision(num_places_);
  double radius, distance;
  int center_index;
  //vmd_.get_params(config, site.type(), &radius, &distance, &center_index);
  ASSERT(config.dimension() == 3, "Euler angles assume 3D.");
  DEBUG("center index " << center_index);
  for (int dim = 0; dim < config.dimension(); ++dim) {
    (*file_.get()) << site.position().coord(dim) << " ";
  }
  const Euler& euler = site.euler();
  (*file_.get()) << euler.phi() << " "
                 << euler.theta() << " "
                 << euler.psi() << std::endl;
}

void FileXYZEuler::write(const std::string file_name,
                    const Configuration& config,
                    const int num_places) const {
  auto file = std::make_shared<std::ofstream>();
  if (append_ == 0) {
    file->open(file_name);
  } else {
    file->open(file_name, std::ofstream::app);
  }
  const Domain& domain = config.domain();
  (*file.get()) << config.group_selects()[group_index_].num_sites() << std::endl
    << "-1 ";
  (*file.get()) << std::setprecision(num_places);
  for (int dim = 0; dim < domain.dimension(); ++dim) {
    (*file.get()) << domain.side_length(dim) << " ";
  }
  (*file.get()) << domain.xy() << " "
    << domain.xz() << " "
    << domain.yz() << " "
    << std::endl;
  PrinterXYZEuler printer(file);
  VisitConfiguration().loop(config, &printer, group_index_);
}

void FileXYZEuler::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2936, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(append_, ostr);
}

FileXYZEuler::FileXYZEuler(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2936, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&append_, istr);
}
}  // namespace feasst
