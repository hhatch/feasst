#include <vector>
#include <string>
#include "system/include/potential_factory.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"

namespace feasst {

void PotentialFactory::add(std::shared_ptr<Potential> potential) {
  potentials_.push_back(potential);
}

void PotentialFactory::set(const int index,
                           std::shared_ptr<Potential> potential) {
  potentials_[index] = potential;
}

void PotentialFactory::precompute(Configuration * config) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    DEBUG("precomputing " << potential->model().class_name());
    potential->precompute(config);
  }
}

void PotentialFactory::precompute(const int index, Configuration * config) {
  potentials_[index]->precompute(config);
}

double PotentialFactory::energy(Configuration * config) {
  double en = 0;
  int index = 0;
  while ((index < num()) && (en < NEAR_INFINITY/10.)) {
    DEBUG("potential index: " << index);
    en += potentials_[index]->energy(config);
    ++index;
  }
  DEBUG("en " << en);
  DEBUG(str());
  return en;
}

double PotentialFactory::energy(const Select& select, Configuration * config) {
  double en = 0;
  int index = 0;
  while ((index < static_cast<int>(potentials_.size())) and
         (en < NEAR_INFINITY)) {
    en += potentials_[index]->energy(select, config);
    ++index;
  }
  DEBUG("en " << en);
  DEBUG(str());
  return en;
}

std::vector<double> PotentialFactory::stored_energy_profile() const {
  std::vector<double> en;
  for (const std::shared_ptr<Potential> potential : potentials_) {
    en.push_back(potential->stored_energy());
  }
  return en;
}

double PotentialFactory::stored_energy() const {
  std::vector<double> en = stored_energy_profile();
  return std::accumulate(en.begin(), en.end(), 0.);
}

std::string PotentialFactory::str() const {
  std::stringstream ss;
  ss << "PotentialFactory: ";
  for (const std::shared_ptr<Potential> potential : potentials_) {
    ss << potential->stored_energy() << " ";
  }
  return ss.str();
}

void PotentialFactory::revert(const Select& select) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->revert(select);
  }
}

void PotentialFactory::finalize(const Select& select) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->finalize(select);
  }
}

void PotentialFactory::serialize(std::ostream& sstr) const {
  feasst_serialize_version(8655, sstr);
  feasst_serialize(potentials_, sstr);
}

PotentialFactory::PotentialFactory(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 8655, "unrecognized verison: " << version);
  feasst_deserialize(&potentials_, sstr);
}

void PotentialFactory::load_cache(const bool load) {
  for (std::shared_ptr<Potential> potential : potentials_) {
    potential->load_cache(load);
  }
}

void PotentialFactory::unload_cache(const PotentialFactory& factory) {
  ASSERT(num() == factory.num(), "size mismatch");
  for (int ip = 0; ip < num(); ++ip) {
    potentials_[ip]->unload_cache(*factory.potentials_[ip]);
  }
}

void PotentialFactory::check() const {
  for (const std::shared_ptr<Potential> potential : potentials_) {
    potential->check();
  }
}

void PotentialFactory::synchronize_(const PotentialFactory& factory,
  const Select& perturbed) {
  for (int index = 0; index < num(); ++index) {
    potentials_[index]->synchronize_(*factory.potentials()[index], perturbed);
  }
}

void PotentialFactory::change_volume(const double delta_volume,
    const int dimension) {
  for (int index = 0; index < num(); ++index) {
    potentials_[index]->change_volume(delta_volume, dimension);
  }
}

}  // namespace feasst
