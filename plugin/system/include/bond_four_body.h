
#ifndef FEASST_CONFIGURATION_BOND_FOUR_BODY_H_
#define FEASST_CONFIGURATION_BOND_FOUR_BODY_H_

#include <map>
#include <string>
#include <memory>
#include "math/include/position.h"
#include "configuration/include/bond.h"

namespace feasst {

/**
  A four body bond is defined by four sites, i - j - k - l.
 */
class BondFourBody {
 public:
  BondFourBody() {}
  virtual double energy(
      const Position& ri,
      const Position& rj,
      const Position& rk,
      const Position& rl,
      const Dihedral& dihedral) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondFourBody> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondFourBody> >& deserialize_map();
  std::shared_ptr<BondFourBody> deserialize(std::istream& istr);
  virtual ~BondFourBody() {}

 protected:
  std::string class_name_ = "BondFourBody";

  void serialize_bond_four_body_(std::ostream& ostr) const;
  explicit BondFourBody(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_FOUR_BODY_H_
