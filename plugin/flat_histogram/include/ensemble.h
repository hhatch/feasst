
#ifndef FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
#define FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_

#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/macrostate.h"

namespace feasst {

class FlatHistogram;

/**
 */
class Ensemble {
 public:
  /// Set the macrostate probability distribution.
  void set_ln_prob(const LnProbability& ln_prob) { ln_prob_ = ln_prob; }

  /// Return the stored LnProbability
  const LnProbability& ln_prob() const { return ln_prob_; }
  
  /// Determine min and max indices for a given phase
  /// Return -1 if no phase boundary.
  void phase_boundary(const LnProbability& ln_prob,
      const int phase, int * min, int * max) const;

  /// Same as above but with the ln_prob_ contained in this class.
  void phase_boundary(const int phase, int * min, int * max) const {
    phase_boundary(ln_prob(), phase, min, max); }

  /// Store and return a reweighted macrostate distribution.
  LnProbability reweight(const LnProbability& ln_prob,
    const Macrostate& macrostate,
    /**
      Change in conjugate variable. For example, if reweighting in number
      of particles using the grand canonical ensemble, the change in the
      thermodynamic conjugate is \f$\Delta(\beta\mu\)f$.
      If reweighting in potential energy in the microcanonical, the change
      in the thermodynamic conjugate is \f$\Delta(-\beta)\f$.
     */
    const double delta_conjugate);

  /// Same as above, but used the stored ln_prob.
  LnProbability reweight(const Macrostate& macrostate,
    const double delta_conjugate) {
    return reweight(ln_prob(), macrostate, delta_conjugate); }

  /// Same as above, but obtain the ln_prob and macrostate from FlatHistogram.
  LnProbability reweight(const FlatHistogram& flat_hist,
    const double delta_conjugate);

  /// Return the ensemble averaged property from a list of properties averaged
  /// at each macrostate.
  double average(const LnProbability& ln_prob,
                 // HWH make macrostate_averages work with bin averages.
                 const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const;

  /// Same as above but using the current ln_prob.
  double average(const std::vector<double>& macrostate_averages,
                 const int phase = 0) const {
    return average(ln_prob(), macrostate_averages, phase); }

  /// Return the average macrostate
  double average(const LnProbability& ln_prob,
      const Macrostate& macrostate,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;

  /// Same as above but using the current ln_prob.
  double average(const Macrostate& macrostate,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const {
    return average(ln_prob(), macrostate, phase);
  }

 private:
  LnProbability ln_prob_;
};

class GrandCanonicalEnsemble : public Ensemble {
 public:
  
  // HWH add reference to Shen/Errington
  /// Return \f$\beta PV\f$, the pressure times inverse temperature and volume.
  double betaPV(const LnProbability& ln_prob,
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;

  /// Same as above but using the current ln_prob.
  double betaPV(const int phase = 0) const {
    return betaPV(ln_prob(), phase); }

};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
