
#ifndef FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Perform a number of attempts.
 */
class Run : public Action {
 public:
  /**
    args:
    - num_attempts: run this many Trial attempts (default: -1. e.g., None)
    - until_num_particles: run until this many particles (default: -1. e.g., None)
   */
  explicit Run(argtype args = argtype());
  explicit Run(argtype * args);
  void perform(MonteCarlo * mc);
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Run>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Run>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Run(std::istream& istr);
  virtual ~Run() {}

 private:
  int num_attempts_;
  int until_num_particles_;
};

inline std::shared_ptr<Run> MakeRun(argtype args = argtype()) {
  return std::make_shared<Run>(args);
}

/**
  Remove a Trial.
 */
class RemoveTrial : public Action {
 public:
  /**
    args:
    - index: index of trial to remove, in order of trials added.
      If -1, do nothing. (default: -1).
   */
  explicit RemoveTrial(argtype args = argtype());
  explicit RemoveTrial(argtype * args);
  void perform(MonteCarlo * mc);
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveTrial>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveTrial>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveTrial(std::istream& istr);
  virtual ~RemoveTrial() {}

 private:
  int index_;
};

inline std::shared_ptr<RemoveTrial> MakeRemoveTrial(argtype args = argtype()) {
  return std::make_shared<RemoveTrial>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
