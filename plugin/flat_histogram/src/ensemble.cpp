#include <cmath>
#include "utils/include/serialize.h"  // deep_copy
#include "flat_histogram/include/ensemble.h"

namespace feasst {

Ensemble::Ensemble(const FlatHistogram& flat_hist) {
  flat_hist_ = deep_copy(flat_hist);//flat_hist;
  ASSERT(ln_prob_original().size() == macrostate().histogram().size(),
    "ln_prob(" << ln_prob_original().size() << ") or macrostate (" << 
    macrostate().histogram().size() << " ) are not of same size.");
  ln_prob_ = ln_prob_original();
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
    ln_prob_.add(macro, macrostate().histogram().center_of_bin(macro)
                 *delta_conjugate);
  }
  ln_prob_.normalize();
  return ln_prob_;
}

double Ensemble::original_conjugate() const { FATAL("not implemented"); }

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
    average += macrostate().value(bin)*std::exp(ln_prob().value(bin));
  }
  return average/ln_prob().sum_probability(min, max);
}

GrandCanonicalEnsemble::GrandCanonicalEnsemble(
    const FlatHistogram& flat_histogram) : Ensemble(flat_histogram) {
}

double GrandCanonicalEnsemble::original_conjugate() const {
  return flat_histogram().beta_mu(); 
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
