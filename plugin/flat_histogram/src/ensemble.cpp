#include <cmath>
#include "utils/include/serialize.h"
#include "flat_histogram/include/ensemble.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"

namespace feasst {

void Ensemble::init_(const Histogram& macrostates,
                     const LnProbability& ln_prob) {
  ln_prob_original_ = ln_prob;
  macrostates_ = macrostates;
  ASSERT(ln_prob_original().size() == macrostates_.size(),
    "ln_prob(" << ln_prob_original().size() << ") or macrostate (" <<
    macrostates_.size() << " ) are not of same size.");
  ln_prob_ = ln_prob_original();
}

Ensemble::Ensemble(const FlatHistogram& flat_hist)
  : Ensemble(flat_hist.macrostate().histogram(),
             flat_hist.ln_prob()) {}

Ensemble::Ensemble(const Clones& clones) {
  Histogram macrostates;
  LnProbability ln_prob = clones.ln_prob(&macrostates);
  init_(macrostates, ln_prob);
}

void Ensemble::phase_boundary(const int phase, int * min, int * max) const {
  std::vector<int> mins = ln_prob().minima();
  const double num_min = static_cast<int>(mins.size());
  if (num_min == 0) {
    *min = 0;
    *max = ln_prob().size() - 1;
  } else if (num_min == 1) {
    if (phase == 0) {
      *min = 0;
      *max = mins[0];
    } else if (phase == 1) {
      *min = mins[0];
      *max = ln_prob().size() - 1;
    } else {
      ERROR("unrecognized phase: " << phase);
    }
  } else {
    ERROR("multiple minima: " << num_min << " not implemented");
  }
}

const LnProbability& Ensemble::reweight(const double delta_conjugate) {
  delta_conjugate_ = delta_conjugate;
  ln_prob_ = ln_prob_original();
  for (int macro = 0; macro < ln_prob_.size(); ++macro) {
    ln_prob_.add(macro, macrostates().center_of_bin(macro)
                 *delta_conjugate);
  }
  ln_prob_.normalize();
  return ln_prob_;
}

double Ensemble::average(const std::vector<double>& macrostate_averages,
     const int phase) const {
  ASSERT(ln_prob().size() == static_cast<int>(macrostate_averages.size()),
    "size mismatch: ln_prob:" << ln_prob().size() <<
    " macro:" << macrostate_averages.size());
  int min, max;
  phase_boundary(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate_averages[bin]*std::exp(ln_prob().value(bin));
  }
  return average/ln_prob().sum_probability(min, max);
}

double Ensemble::average_macrostate(const int phase) const {
  int min, max;
  phase_boundary(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostates().center_of_bin(bin)*std::exp(ln_prob().value(bin));
  }
  return average/ln_prob().sum_probability(min, max);
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(const Histogram& macrostates,
  const LnProbability& ln_prob,
  const double conjugate) : Ensemble(macrostates, ln_prob) {
  original_conjugate_ = conjugate;
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(
    const FlatHistogram& flat_histogram) : Ensemble(flat_histogram) {
  original_conjugate_ = flat_histogram.beta_mu();
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(
    const Clones& clones) : Ensemble(clones) {
  original_conjugate_ = clones.clone(0).criteria().beta_mu();
}

double GrandCanonicalEnsemble::betaPV(const int phase) const {
  int min, max;
  phase_boundary(phase, &min, &max);
  return -ln_prob().value(0) + std::log(ln_prob().sum_probability(min, max));
}

//double beta() const;
//  { return beta_; }
//  double beta_mu() const { return beta_chemical potential

}  // namespace feasst
