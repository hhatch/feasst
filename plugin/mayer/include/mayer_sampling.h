
#ifndef FEASST_MAYER_MAYER_SAMPLING_H_
#define FEASST_MAYER_MAYER_SAMPLING_H_

#include "monte_carlo/include/criteria.h"
#include "math/include/accumulator.h"

namespace feasst {

class Random;

/**
  Mayer-sampling Monte Carlo acceptance criteria (see
  https://doi.org/10.1103/PhysRevLett.92.220601).
 */
class MayerSampling : public Criteria {
 public:
  MayerSampling();

  bool is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) override;

  /// Return the mayer ensemble of the full potential.
  const Accumulator& mayer() const { return mayer_; }

  /// Return the mayer ensemble of the reference potential.
  const Accumulator& mayer_ref() const { return mayer_ref_; }

  /// Return the ratio of the second virial coefficient of the full potential
  /// to the second virial coefficient of the reference potential.
  double second_virial_ratio() const;

  /**
    Return the block standard deviation of the second_virial_ratio.
    This is computed using the error propogation formula variance formula:

    \f$f = mayer/mayer_ref\f$

    \f$\sigma_f = \sqrt{(\sigma_{mayer}/mayer_ref)^2 + (f\sigma_{mayer_ref}/mayer_ref)^2}\f$
   */
  double second_virial_ratio_block_stdev() const;

  std::shared_ptr<Criteria> create(std::istream& istr) const override {
    return std::make_shared<MayerSampling>(istr); }
  std::shared_ptr<Criteria> create(argtype * args) const override {
    return std::make_shared<MayerSampling>(); }
  void serialize(std::ostream& ostr) const override;
  explicit MayerSampling(std::istream& istr);
//  explicit MayerSampling(const Criteria& criteria);
  ~MayerSampling() {}

 private:
  double f12old_ = -1.;
  double f12ref_ = -1.;
  Accumulator mayer_;
  Accumulator mayer_ref_;
};

inline std::shared_ptr<MayerSampling> MakeMayerSampling() {
  return std::make_shared<MayerSampling>();
}

}  // namespace feasst

#endif  // FEASST_MAYER_MAYER_SAMPLING_H_
