
#ifndef FEASST_PATCH_FILE_XYZ_PATCH_H_
#define FEASST_PATCH_FILE_XYZ_PATCH_H_

#include <string>
#include <fstream>
#include <sstream>
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"

namespace feasst {

// Utility class to print XYZ files from selection.
class PrinterXYZEuler : public LoopConfigOneBody {
 public:
  PrinterXYZEuler(std::shared_ptr<std::ofstream> file, const int num_places = 8);
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override;

 private:
  int num_places_;
  std::shared_ptr<std::ofstream> file_;
};

/**
  Similar to FileXYZ, except include Euler angles also.
 */
class FileXYZEuler {
 public:
  /**
    args:
    - group_index: print the coordinates of this group index only (default: 0).
    - append: append file output if set to true.
      Do not append if false (default: "false").
   */
  explicit FileXYZEuler(argtype args = argtype());
  explicit FileXYZEuler(argtype * args);

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

  bool append() const { return append_; }
  void serialize(std::ostream& ostr) const;
  explicit FileXYZEuler(std::istream& istr);

 private:
  int group_index_;
  bool append_;
};

inline std::shared_ptr<FileXYZEuler> MakeFileXYZEuler(argtype args = argtype()) {
  return std::make_shared<FileXYZEuler>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_FILE_XYZ_PATCH_H_
