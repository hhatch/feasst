#include "monte_carlo/include/criteria.h"
#include "steppers/include/criteria_updater.h"

namespace feasst {

void CriteriaUpdater::update(Criteria * criteria,
  System * system,
  Random * random,
  TrialFactory * trial_factory) { criteria->update(); }

} // namespace feasst
