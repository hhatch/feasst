
#ifndef FEASST_CONFIGURATION_BOND_THREE_BODY_H_
#define FEASST_CONFIGURATION_BOND_THREE_BODY_H_

#include <map>
#include <string>
#include <memory>
#include "math/include/position.h"
#include "configuration/include/bond.h"

namespace feasst {

/**
  A three body bond is defined by three sites, 0 - 1 - 2.
  The relative vector r01 = r0 - r1 points from r1 to r0.
  The relative vector r21 = r2 - r1 points from r1 to r2.
 */
class BondThreeBody {
 public:
  BondThreeBody() {}
  virtual double energy(const Position& relative01, const Position& relative21,
    const Bond& angle) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondThreeBody> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondThreeBody> >& deserialize_map();
  std::shared_ptr<BondThreeBody> deserialize(std::istream& istr);
  virtual ~BondThreeBody() {}

 protected:
  std::string class_name_ = "BondThreeBody";

  void serialize_bond_three_body_(std::ostream& ostr) const;
  explicit BondThreeBody(std::istream& istr);
};

class AngleModel : public BondThreeBody {
 public:
  AngleModel() {}
  double energy(const Position& relative01, const Position& relative21,
    const Bond& angle) const override;
  virtual double energy(const double theta, const Bond& angle) const = 0;
  explicit AngleModel(std::istream& istr) : BondThreeBody(istr) {}
  virtual ~AngleModel() {}
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_THREE_BODY_H_
