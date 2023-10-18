
#ifndef FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_
#define FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_

#include <string>
#include <sstream>
#include <memory>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/collection_matrix.h"
#include "steppers/include/seek_analyze.h"

namespace feasst {

class Checkpoint;
class Histogram;

/**
  Container for holding a group of FlatHistogram MonteCarlo simulations
  where each of the collection matricies can be spliced together.

  Beware attempting to splice windows with different move sets (trials/weights).
  This can cause issues in the calculation of ln_prob.

  Adjust bounds based on the number of iterations.
  Thus, iteratively give a window with less iterations macrostates from a
  window with more iterations, while possible.
  Left and right most windows can completely abandon completed macrostates.
  Also, left and right most don't lose macrostates if the entire window
  has already reached its completion criteria.
 */
class CollectionMatrixSplice {
 public:
  //@{
  /** @name Arguments
    - min_window_size: minimum size of window boundaries during adjustment.
      If -1, do not adjust window boundaries (default: 5).
    - hours_per: hours per bounds adjustment, checkpoint and writing the
      combined ln_prob (default: 0.01).
    - ln_prob_file: file name for the combined ln_prob, if not empty (default: empty).
    - ln_prob_file_append: if true, append to ln_prob_file (default: false).
    - num_adjust_per_write: number of adjustments per writing of ln_prob (default: 1)
    - bounds_file: file name for periodic output of bounds, if not empty
      (default: empty).
   */
  explicit CollectionMatrixSplice(argtype args = argtype());
  explicit CollectionMatrixSplice(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Add a MonteCarlo.
  void add(std::shared_ptr<MonteCarlo> mc) { clones_.push_back(mc); }

  /// Set size
  void set_size(const int size) { clones_.resize(size); }

  /// Set a MonteCarlo.
  void set(const int index, std::shared_ptr<MonteCarlo> mc) {
    clones_[index] = mc; }

//  /// Create all clones via deep copy of given object and given windows.
//  CollectionMatrixSplice(const MonteCarlo& mc, const Window& window);

  /// Return the number of clones.
  int num() const { return static_cast<int>(clones_.size()); }

  /// Return a read-only clone.
  const MonteCarlo& clone(const int index) const;

  /// Return a writable clone.
  MonteCarlo * get_clone(const int index);

  /// Return the clones.
  std::vector<std::shared_ptr<MonteCarlo> > get_clones() { return clones_; }

  /// Set the Checkpoint
  void set(std::shared_ptr<Checkpoint> checkpoint);

  /// Return the FlatHistogram of a given clone index.
  const FlatHistogram& flat_histogram(const int index) const;

  /// Return the CollectionMatrix of a given clone index.
  const CollectionMatrix& collection_matrix(const int index) const;

  /// Loop through all clones and swap bounds based on number of interations.
  void adjust_bounds();

  /// Write the bounds to file, if not empty.
  void write_bounds(
    /// optionally, print the header.
    const bool header = false);

  /// Return true if all clones are complete
  bool are_all_complete() const;

  /// Run for a number of hours.
  void run(const double hours);

  /// Run until all are complete.
  void run_until_all_are_complete();

  /// Return the complete collection matrix.
  CollectionMatrix collection_matrix() const;

  /// Return the complete probability distribution.
  LnProbability ln_prob() const;

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit CollectionMatrixSplice(std::istream& istr);

  /// write to file
  void write(const std::string& file_name) const;

  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  CollectionMatrixSplice deserialize(const std::string str) {
    std::stringstream ss(str);
    return CollectionMatrixSplice(ss);
  }

  //@}
 private:
  std::vector<std::shared_ptr<MonteCarlo> > clones_;
  std::shared_ptr<Checkpoint> checkpoint_;
  int min_window_size_;
  double hours_per_;
  std::string ln_prob_file_;
  bool ln_prob_file_append_;
  int num_adjust_per_write_;
  int num_adjust_since_write_ = 0;
  std::string bounds_file_;
};

/// Construct CollectionMatrixSplice
inline std::shared_ptr<CollectionMatrixSplice> MakeCollectionMatrixSplice(
  argtype args = argtype()) {
  return std::make_shared<CollectionMatrixSplice>(args); }

/// Construct MonteCarlo from Checkpoint file.
std::shared_ptr<CollectionMatrixSplice> MakeCollectionMatrixSplice(
  const std::string file_name);

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_
