
#ifndef FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
 */
class Run : public Action {
 public:
  /**
    args:
    - num_attempts: run this many Trial attempts (default: -1. e.g., None)
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
};

inline std::shared_ptr<Run> MakeRun(argtype args = argtype()) {
  return std::make_shared<Run>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
