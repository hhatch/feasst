
#ifndef FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
#define FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_

#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {

//class FlatHistogram;

/**
  Perform reweighting and ensemble averages using macrostate distributions
 */
class Ensemble {
 public:
  /// Store the original macrostate distribution
  Ensemble(const FlatHistogram& flat_hist);

  /// Return the stored macrostate distribution from the original simulation.
  const LnProbability& ln_prob_original() const { return flat_hist_.ln_prob(); }

  /// Return the LnProbability (that may have been reweighted).
  const LnProbability& ln_prob() const { return ln_prob_; }

  /// Return the stored macrostate from the original simulation.
  const Macrostate& macrostate() const { return flat_hist_.macrostate(); }

  /// Return the stored FlatHistogram from the original simulation.
  const FlatHistogram& flat_histogram() const { return flat_hist_; }

  /// Determine min and max macrostate indices for a given phase
  /// Return -1 if no phase boundary.
  void phase_boundary(const int phase, int * min, int * max) const;

  /// Return the conjugate thermodynamic variable of the macrostate from the
  /// original simulation.
  virtual double original_conjugate() const;

  /**
    Store and return a reweighted macrostate distribution due to a change
    in the conjugate variable.
    For example, if reweighting in number
    of particles using the grand canonical ensemble, the change in the
    thermodynamic conjugate is \f$\Delta(\beta\mu\)f$.
    If reweighting in potential energy in the microcanonical, the change
    in the thermodynamic conjugate is \f$\Delta(-\beta)\f$.
   */
  const LnProbability& reweight(const double delta_conjugate);

  /// Return the conjugate variable corresponding to the current (reweighted)
  /// macrostate distribution.
  double conjugate() const { return original_conjugate() + delta_conjugate_; }

  /// Return the ensemble averaged property from a list of properties averaged
  /// at each macrostate.
  double average(const std::vector<double>& macrostate_averages,
                 /// Select phase by order of macrostate.
                 /// Assumes default method of dividing phase boundary.
                 const int phase = 0) const;

  /// Return the average macrostate
  double average_macrostate(
    /// Select phase by order of macrostate.
    /// Assumes default method of dividing phase boundary.
    const int phase = 0) const;

 private:
  FlatHistogram flat_hist_;
  LnProbability ln_prob_;
  double delta_conjugate_ = 0.; // record reweight change in conjugate
};

/**
  Grand canonical ensemble currently implemented for single component systems.
 */
class GrandCanonicalEnsemble : public Ensemble {
 public:
  /// Store \f$\beta\f$, \f$\mu\f$ and the original ln_prob_.
  GrandCanonicalEnsemble(const FlatHistogram& flat_histogram);
  double original_conjugate() const override;

  /// Return the conjugate to the macrostate, \f$\beta\mu\f$.
  double beta_mu() const { return conjugate(); }

  /// Return \f$\beta PV\f$, the pressure times volume and inverse temperature.
  /// See https://doi.org/10.1063/1.2064628
  double betaPV(
      /// Select phase by order of macrostate.
      /// Assumes default method of dividing phase boundary.
      const int phase = 0) const;

 private:
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_ENSEMBLE_H_
