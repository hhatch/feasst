#include <cmath>
#include "utils/include/debug.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

double Acceptance::ln_metropolis_prob() const {
  ASSERT(!std::isinf(ln_metropolis_prob_), "ln_metropolis_prob_ is inf");
  return ln_metropolis_prob_;
}

void Acceptance::reset() {
  set_ln_metropolis_prob();
  set_reject();
  set_endpoint();
  energy_new_ = 0.;
  energy_old_ = 0.;
  configuration_index_ = 0.;
  std::fill(energy_profile_new_.begin(), energy_profile_new_.end(), 0.);
  std::fill(energy_profile_old_.begin(), energy_profile_old_.end(), 0.);
  macrostate_shift_ = 0;
  macrostate_shift_type_ = 0;
  perturbed_.clear();
  perturbed_.resize(2); // maximum number of configs
}

void Acceptance::add_to_perturbed(const Select& select, const int config) {
  if (config >= static_cast<int>(perturbed_.size())) {
    perturbed_.resize(config + 1);
  }
  perturbed_[config].add(select);
}

void Acceptance::set_perturbed_state(const int state, const int config) {
  if (config >= static_cast<int>(perturbed_.size())) {
    perturbed_.resize(config + 1);
  }
  DEBUG("state " << state);
  perturbed_[config].set_trial_state(state);
}

const Select& Acceptance::perturbed(const int config) const {
  ASSERT(config < static_cast<int>(perturbed_.size()),
    "config: " << config << " >= size:" << perturbed_.size() <<
    "Consider increasing max num config in reset()");
  return perturbed_[config];
}

void Acceptance::add_to_energy_profile_new(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_new_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_new_.resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[i] += energy[i];
  }
}

void Acceptance::subtract_from_energy_profile_new(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_new_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_new_.resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[i] -= energy[i];
  }
}

void Acceptance::add_to_energy_profile_old(const std::vector<double>& energy) {
  const int current_size = static_cast<int>(energy_profile_old_.size());
  const int new_size = static_cast<int>(energy.size());
  if (current_size < new_size) {
    energy_profile_old_.resize(new_size);
  } else {
    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_old_[i] += energy[i];
  }
}

void Acceptance::set_reject(const bool reject) {
  reject_ = reject;
}

}  // namespace feasst
