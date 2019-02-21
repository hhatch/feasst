
#ifndef FEASST_CORE_BOND_H_
#define FEASST_CORE_BOND_H_

#include <vector>
#include "core/include/typed_entity.h"
#include "core/include/properties.h"

namespace feasst {

/**
  Sites within the same particle may be bonded.
  The indices of the sites which are bonded are stored here.
  The type of the bond is used to determine the bond model.
 */
class Bond : public PropertiedEntity, public TypedEntity {
 public:
  /// Return the indices of the sites involved in the bond within a particle.
  std::vector<int> site_indices() const { return site_indicies_; }

  /// Return the indices of the sites involved in the bond within a particle.
  int site(const int index) const { return site_indicies_[index]; }

  /// Return the number of sites in bond.
  int num_sites() const { return static_cast<int>(site_indicies_.size()); }

  /// Add site index.
  void add_site_index(const int index) { site_indicies_.push_back(index); }

 private:
  std::vector<int> site_indicies_;
};

class Angle : public Bond {};
class Dihedral : public Bond {};
class Improper : public Bond {};

// HWH add dihedrals and impropers here

}  // namespace feasst

#endif  // FEASST_CORE_BOND_H_
