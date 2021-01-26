
#ifndef FEASST_CONFIGURATION_FILE_XYZ_H_
#define FEASST_CONFIGURATION_FILE_XYZ_H_

#include <string>
#include <fstream>
#include <sstream>
#include "configuration/include/configuration.h"

namespace feasst {

class FileVMD {
 public:
  FileVMD() {}

  void write(const std::string file_name,
      const Configuration& config,
      const std::string traj_file_name);

  void serialize(std::ostream& ostr) const;
  explicit FileVMD(std::istream& istr);
};

// Note HWH: best not to assume this second line format by default
/**
  The XYZ file format has no formal standard: https://en.wikipedia.org/wiki/XYZ_file_format
  It is important to remember that FEASST as its own variant.

  The first line is the number of sites, n.
  The second line is of the format [id lx ly lz xy xz yz], where id is a
  placeholder for order parameters or macrostates,
  lx is the box length in the x direction,
  ly is the box length in the y direction,
  lz is the box length in the z direction,
  xy is the xy domain tilt (see Domain),
  xz is the xz domain tilt, and yz is the yz domain tilt.
  The following n lines are in the format [id x y z],
  where id is the unique site type and x, y, z are the Cartesian coordinates.
 */
class FileXYZ {
 public:
  FileXYZ() {
    set_group_index();
    set_append();
  }

  /**
    Load the xyz file with file_name into the configuration.
    Note that this function does not read the domain tilts xy, xz and yz.
    This function also does not read the site types.
    Thus, particles should be added to the system in the desired order.
    If no particles in the configuration, use the first particle type.
   */
  void load(const std::string file_name, Configuration * config) const;

  /// Write the configuration to file_name in xyz format.
  /// If the simulation is 2D, simply writes z as 0.
  void write(const std::string file_name,
             const Configuration& config,
             /// Number of decimal places
             const int num_decimal_places = 8) const;

  void write_for_vmd(const std::string file_name,
             const Configuration& config) const;

  void set_group_index(const int index = 0) { group_index_ = index; }

  /// By default, do not append
  void set_append(const int append = 0) { append_ = append; }

  void serialize(std::ostream& ostr) const;
  explicit FileXYZ(std::istream& istr);

 private:
  int group_index_;
  int append_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_XYZ_H_
