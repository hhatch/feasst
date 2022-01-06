#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "monte_carlo/include/trial_factory.h"
#include "math/include/random.h"

namespace feasst {

TrialFactory::TrialFactory(argtype args) : TrialFactory(&args) {
  check_all_used(args);
}
TrialFactory::TrialFactory(argtype * args) : Trial(args) {
  class_name_ = "TrialFactory";
}

class MapTrialFactory {
 public:
  MapTrialFactory() {
    auto obj = MakeTrialFactory();
    obj->deserialize_map()["TrialFactory"] = obj;
  }
};

static MapTrialFactory mapper_ = MapTrialFactory();

void TrialFactory::add(std::shared_ptr<Trial> trial) {
  trials_.push_back(trial);
  update_cumul_prob_();
}

void TrialFactory::update_cumul_prob_() {
  std::vector<double> weights;
  for (std::shared_ptr<Trial> trial : trials_) {
    weights.push_back(trial->weight());
  }
  if (weights.size() > 0) {
    cumulative_probability_ = feasst::cumulative_probability(weights);
  }
  //std::stringstream ss;
  //ss << trials_.back()->class_name()"trial" << num() - 1;
  // timer_.add(trials_.back()->class_name());
}

void TrialFactory::remove(const int index) {
  trials_.erase(trials_.begin() + index);
  update_cumul_prob_();
}

int TrialFactory::random_index(Random * random) {
  ASSERT(num() > 0, "no trials to select");
  return random->index_from_cumulative_probability(cumulative_probability_);
}

bool TrialFactory::attempt(
    Criteria* criteria,
    System * system,
    int trial_index,
    Random * random) {
  increment_num_attempts();
  // timer_.start(0);
  if (num() == 0) return false;
  if (trial_index == -1) {
    trial_index = random_index(random);
  }
  last_index_ = trial_index;
  //timer_.start(index + 1);  // +1 for "other"
  const bool accepted = trials_[trial_index]->attempt(criteria, system, random);
  //timer_.end();
  if (accepted) increment_num_success_();
  return accepted;
}

void TrialFactory::revert(const int index,
                          const bool accepted,
                          const bool auto_rejected,
                          System * system) {
  trials_[index]->revert(index, accepted, auto_rejected, system);
  if (accepted) {
    decrement_num_success_();
  }
  decrement_num_attempts_();
}

void TrialFactory::imitate_trial_rejection_(const int index,
    const bool auto_reject) {
  DEBUG("index " << index << " " << trials_.size());
  trials_[index]->increment_num_attempts();
  if (auto_reject) {
    trials_[index]->increment_num_auto_reject();
  }
  increment_num_attempts();
}

std::string TrialFactory::status_header() const {
  std::stringstream ss;
  for (int trial = 0; trial < num(); ++trial) {
    ss << trials_[trial]->status_header();
  }
  return ss.str();
}

std::string TrialFactory::status() const {
  std::stringstream ss;
  for (int trial = 0; trial < num(); ++trial) {
    ss << trials_[trial]->status();
  }
  return ss.str();
}

void TrialFactory::delay_finalize() {
  for (int trial = 0; trial < num(); ++trial) {
    trials_[trial]->set_finalize_delayed(true);
  }
}

void TrialFactory::reset_stats() {
  Trial::reset_stats();
  for (int trial = 0; trial < num(); ++trial) {
    trials_[trial]->reset_stats();
  }
}

void TrialFactory::tune() {
  for (int trial = 0; trial < num(); ++trial) {
    trials_[trial]->tune();
  }
}

void TrialFactory::precompute(Criteria * criteria, System * system) {
  for (std::shared_ptr<Trial> trial : trials_) {
    trial->precompute(criteria, system);
  }
}

bool TrialFactory::is_equal(const TrialFactory& factory) const {
  if (num() != factory.num()) {
    INFO("unequal number of trials: " << num() << " "
      << factory.num());
    return false;
  }
  if (!Trial::is_equal(factory)) {
    INFO("here");
    return false;
  }
  for (int it = 0; it < num(); ++it) {
    if (!trials_[it]->is_equal(*(factory.trials_[it]))) {
      INFO("unequal trial" << it);
      return false;
    }
  }
  return true;
}

std::shared_ptr<Trial> TrialFactory::create(std::istream& istr) const {
  return std::make_shared<TrialFactory>(istr);
}

TrialFactory::TrialFactory(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialFactory", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(189 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize(&trials_, istr);
  { int dim1;
    istr >> dim1;
    trials_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstobj((trials_)[index], istr);
      int existing;
      istr >> existing;
      if (existing != 0) {
        trials_[index] = trials_[index]->deserialize(istr);
      }
    }
  }
  feasst_deserialize(&cumulative_probability_, istr);
}

void TrialFactory::serialize_trial_factory_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(189, ostr);
  feasst_serialize(trials_, ostr);
  feasst_serialize(cumulative_probability_, ostr);
}

void TrialFactory::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_factory_(ostr);
}

void TrialFactory::synchronize_(const Trial& trial) {
  Trial::synchronize_(trial);
  for (int itrial = 0; itrial < num(); ++itrial) {
    trials_[itrial]->synchronize_(trial.trial(itrial));
  }
}

std::map<std::string, std::shared_ptr<TrialFactoryNamed> >& TrialFactoryNamed::deserialize_map() {
  static std::map<std::string, std::shared_ptr<TrialFactoryNamed> >* ans =
     new std::map<std::string, std::shared_ptr<TrialFactoryNamed> >();
  return *ans;
}

//void TrialFactoryNamed::serialize(std::ostream& ostr) const {
//  ostr << class_name() << " ";
//  serialize_trial_(ostr);
//}

//std::shared_ptr<TrialFactoryNamed> TrialFactoryNamed::create(std::istream& istr) const {
//  return std::make_shared<TrialFactoryNamed>(istr);
//}

std::shared_ptr<TrialFactoryNamed> TrialFactoryNamed::create(argtype * args) const {
  FATAL("not implemented");
}

//std::shared_ptr<TrialFactoryNamed> TrialFactoryNamed::deserialize(std::istream& istr) {
//  return template_deserialize(deserialize_map(), istr,
//    // true argument denotes rewinding to reread class name
//    // this allows derived class constructor to read class name.
//    true);
//}

std::shared_ptr<TrialFactoryNamed> TrialFactoryNamed::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

}  // namespace feasst
