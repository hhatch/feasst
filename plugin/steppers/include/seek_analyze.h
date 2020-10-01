
#ifndef FEASST_STEPPERS_SEEK_ANALYZE_H_
#define FEASST_STEPPERS_SEEK_ANALYZE_H_

#include <string>
#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

class MonteCarlo;

class AnalyzeData {
 public:
  virtual double get(const Analyze& analyze) const = 0;
};

/// Obtain Accumulator average in Analyze.
class AccumulatorAverage : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override {
    return analyze.accumulator().average();
  }
};

/// Obtain Accumulator sum in Analyze.
class AccumulatorSum : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override {
    return analyze.accumulator().sum();
  }
};

/// Obtain Accumulator sum of squared in Analyze.
class AccumulatorSumOfSquared : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override {
    return analyze.accumulator().sum_of_squared();
  }
};

/**
  Find Analyze with class name.
 */
class SeekAnalyze {
 public:
  /**
    Return the indices, where the first is mc.analyze index.
    If inside AnalyzeFactory, the second is the index of the factory.
    Otherwise, the second index is -1.
    Only the first match is returned.
   */
  std::vector<int> index(const std::string class_name,
                         const MonteCarlo& mc) const;

  /**
    For multistate Analyze with given class_name,
    Return average Accumulator as a function of state.
   */
  std::vector<double> multistate_average(
    const std::string class_name,
    const MonteCarlo& mc,
    /// Optionally specify where to get data from Analyze
    const AnalyzeData& get = AccumulatorAverage()) const;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SEEK_ANALYZE_H_
