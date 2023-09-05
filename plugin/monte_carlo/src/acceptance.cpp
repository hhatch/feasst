#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/utils.h"
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
  energy_new_.resize(1);
  energy_new_[0] = 0.;
  energy_old_.resize(1);
  energy_old_[0] = 0.;
  energy_ref_.resize(1);
  energy_ref_[0] = 0.;
  configuration_indices_.resize(1);
  configuration_indices_[0] = 0;
  resize(1, 0, &energy_profile_new_);
  fill(0., &energy_profile_new_);
  resize(1, 0, &energy_profile_old_);
  fill(0., &energy_profile_old_);
  macrostate_shift_.resize(1);
  macrostate_shift_[0] = 0;
  macrostate_shift_type_.resize(1);
  macrostate_shift_type_[0] = 0.;
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

const std::vector<double>& Acceptance::energy_profile_new(const int config) const {
  ASSERT(config < static_cast<int>(energy_profile_new_.size()),
    "config:" << config << " >= size:" << energy_profile_new_.size());
  return energy_profile_new_[config];
}

void Acceptance::set_energy_profile_new(const std::vector<double>& energy, const int config) {
  if (config == static_cast<int>(energy_profile_new_.size())) {
    energy_profile_new_.resize(energy_profile_new_.size()+1);
  }
  energy_profile_new_[config] = energy;
}

void Acceptance::add_to_energy_profile_new(const std::vector<double>& energy,
    const int config) {
  DEBUG("config " << config);
  if (config == static_cast<int>(energy_profile_new_.size())) {
    energy_profile_new_.resize(energy_profile_new_.size()+1);
  }
  const int current_size = static_cast<int>(energy_profile_new_[config].size());
  DEBUG("current_size " << current_size);
  const int new_size = static_cast<int>(energy.size());
  DEBUG("new_size " << new_size);
  if (current_size < new_size) {
    energy_profile_new_[config].resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[config][i] += energy[i];
  }
}

void Acceptance::subtract_from_energy_profile_new(const std::vector<double>& energy,
    const int config) {
  DEBUG("config " << config);
  const int current_size = static_cast<int>(energy_profile_new_[config].size());
  DEBUG("current_size " << current_size);
  const int new_size = static_cast<int>(energy.size());
  DEBUG("new_size " << new_size);
  if (current_size < new_size) {
    energy_profile_new_[config].resize(new_size);
//  } else {
//    ASSERT(current_size == new_size, "err");
  }
  for (int i = 0; i < new_size; ++i) {
    energy_profile_new_[config][i] -= energy[i];
  }
}

const std::vector<double>& Acceptance::energy_profile_old(const int config) const {
  ASSERT(config < static_cast<int>(energy_profile_old_.size()),
    "config:" << config << " >= size:" << energy_profile_old_.size());
  return energy_profile_old_[config];
}

void Acceptance::set_energy_profile_old(const std::vector<double>& energy, const int config) {
  energy_profile_old_[config] = energy;
}

void Acceptance::add_to_energy_profile_old(const std::vector<double>& energy,
    const int config) {
  if (config == static_cast<int>(energy_profile_old_.size())) {
    energy_profile_old_.resize(energy_profile_old_.size()+1);
  }
  const int current_size = static_cast<int>(energy_profile_old_[config].size());
  const int old_size = static_cast<int>(energy.size());
  if (current_size < old_size) {
    energy_profile_old_[config].resize(old_size);
//  } else {
//    ASSERT(current_size == old_size, "err");
  }
  for (int i = 0; i < old_size; ++i) {
    energy_profile_old_[config][i] += energy[i];
  }
}

void Acceptance::set_reject(const bool reject) {
  reject_ = reject;
}

double Acceptance::energy_new(const int config) const {
  return energy_new_[config];
}

void Acceptance::set_energy_new(const double energy, const int config) {
  energy_new_[config] = energy;
}

void Acceptance::add_to_energy_new(const double energy, const int config) {
  energy_new_[config] += energy;
}

double Acceptance::energy_old(const int config) const {
  return energy_old_[config];
}

void Acceptance::set_energy_old(const double energy, const int config) {
  energy_old_[config] = energy;
}

void Acceptance::add_to_energy_old(const double energy, const int config) {
  energy_old_[config] += energy;
}

}  // namespace feasst
