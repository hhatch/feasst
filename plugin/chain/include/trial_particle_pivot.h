#ifndef FEASST_CHAIN_TRIAL_PARTICLE_PIVOT_H_
#define FEASST_CHAIN_TRIAL_PARTICLE_PIVOT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Rigidly pivot a random particle about one of its sites.
class TrialParticlePivot : public TrialMove {
 public:
  explicit TrialParticlePivot(argtype args = argtype());
  explicit TrialParticlePivot(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialParticlePivot>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialParticlePivot>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialParticlePivot(std::istream& istr);
  virtual ~TrialParticlePivot() {}
};

inline std::shared_ptr<TrialParticlePivot> MakeTrialParticlePivot(
  argtype args = argtype()) {
  return std::make_shared<TrialParticlePivot>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PARTICLE_PIVOT_H_
