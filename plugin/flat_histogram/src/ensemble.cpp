#include <cmath>
#include "utils/include/serialize.h"  // deep_copy
#include "flat_histogram/include/ensemble.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {


LnProbability Ensemble::reweight(const LnProbability& ln_prob,
    const Macrostate& macrostate,
    const double delta_conjugate) {
  LnProbability lnpirw = deep_copy(ln_prob);
  for (int macro = 0; macro < lnpirw.size(); ++macro) {
    lnpirw.add(macro, macrostate.histogram().center_of_bin(macro)
               *delta_conjugate);
  }
  lnpirw.normalize();
  ln_prob_ = lnpirw;
  return lnpirw;
}

LnProbability Ensemble::reweight(const FlatHistogram& flat_hist,
    const double delta_conjugate) {
  return reweight(flat_hist.ln_prob(), flat_hist.macrostate(), delta_conjugate);
}

void Ensemble::phase_boundary(const LnProbability& ln_prob,
    const int phase, int * min, int * max) const {
  std::vector<int> mins = ln_prob.minima();
  const double num_min = static_cast<int>(mins.size());
  if (num_min == 0) {
    *min = 0;
    *max = ln_prob.size() - 1;
  } else if (num_min == 1) {
    if (phase == 0) {
      *min = 0;
      *max = mins[0];
    } else if (phase == 1) {
      *min = mins[0];
      *max = ln_prob.size() - 1;
    } else {
      ERROR("unrecognized phase: " << phase);
    }
  } else {
    ERROR("multiple minima: " << num_min << " not implemented");
  }
}

double GrandCanonicalEnsemble::betaPV(const LnProbability& ln_prob,
    const int phase) const {
  int min, max;
  phase_boundary(ln_prob, phase, &min, &max);
  return -ln_prob.value(0) + std::log(ln_prob.sum_probability(min, max));
}

double Ensemble::average(const LnProbability& ln_prob,
     const std::vector<double>& macrostate_averages,
     const int phase) const {
  ASSERT(ln_prob.size() == static_cast<int>(macrostate_averages.size()),
    "size mismatch: ln_prob:" << ln_prob.size() <<
    " macro:" << macrostate_averages.size());
  int min, max;
  phase_boundary(ln_prob, phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate_averages[bin]*std::exp(ln_prob.value(bin));
  }
  return average/ln_prob.sum_probability(min, max);
}

double Ensemble::average(const LnProbability& ln_prob,
    const Macrostate& macrostate,
    const int phase) const {
  int min, max;
  phase_boundary(phase, &min, &max);
  double average = 0.;
  for (int bin = min; bin < max + 1; ++bin) {
    average += macrostate.value(bin)*std::exp(ln_prob.value(bin));
  }
  return average/ln_prob.sum_probability(min, max);
}

}  // namespace feasst
