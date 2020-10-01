
#ifndef FEASST_STEPPERS_SEEK_ANALYZE_H_
#define FEASST_STEPPERS_SEEK_ANALYZE_H_

#include <string>
#include <vector>

namespace feasst {

class MonteCarlo;

/**
  Find Analyze with class name.
 */
class SeekAnalyze {
 public:
  /// Return the indices, where the first is mc.analyze index.
  /// If inside AnalyzeFactory, the second is the index of the factory.
  /// Otherwise, the second index is -1.
  /// Only the first match is returned.
  std::vector<int> index(const std::string class_name,
                         const MonteCarlo& mc) const;

  /// For multistate Analyze with given class_name,
  /// Return average Accumulator as a function of state.
  std::vector<double> multistate_average(const std::string class_name,
                                         const MonteCarlo& mc) const;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SEEK_ANALYZE_H_
