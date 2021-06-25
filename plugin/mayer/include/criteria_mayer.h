
#ifndef FEASST_MAYER_CRITERIA_MAYER_H_
#define FEASST_MAYER_CRITERIA_MAYER_H_

#include "monte_carlo/include/criteria.h"
#include "math/include/accumulator.h"

namespace feasst {

class Random;

// HWH Rename to MayerSampling
/**
  Mayer-sampling Monte Carlo acceptance criteria.
  Currently assumes hard sphere reference.
  HWH add reference.
 */
class CriteriaMayer : public Criteria {
 public:
  CriteriaMayer() : Criteria() {}

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override;

  double second_virial() const;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<CriteriaMayer>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<CriteriaMayer>(); }
  void serialize(std::ostream& ostr) const override;
  CriteriaMayer(std::istream& istr);
  ~CriteriaMayer() {}

 private:
  const std::string class_name_ = "CriteriaMayer";
  double f12old_ = -1.;
  double f12ref_ = -1.;
  Accumulator mayer_;
  Accumulator mayer_ref_;
};

inline std::shared_ptr<CriteriaMayer> MakeCriteriaMayer() {
  return std::make_shared<CriteriaMayer>();
}

}  // namespace feasst

#endif  // FEASST_MAYER_CRITERIA_MAYER_H_
