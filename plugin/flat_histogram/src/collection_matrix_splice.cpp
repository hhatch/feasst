#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
//#ifdef _WIN32
//  #include <Windows.h>  // sleep
//#else
//  #include <unistd.h>  // sleep
//#endif
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <fstream>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/checkpoint.h"
#include "math/include/constants.h"
#include "math/include/histogram.h"
#include "math/include/utils_math.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/collection_matrix_splice.h"
#include "monte_carlo/include/run.h"

namespace feasst {

CollectionMatrixSplice::CollectionMatrixSplice(argtype * args) {
  min_window_size_ = integer("min_window_size", args, 2);
  hours_per_ = dble("hours_per", args, 0.01);
  ln_prob_file_ = str("ln_prob_file", args, "");
  ln_prob_file_append_ = boolean("ln_prob_file_append", args, false);
}
CollectionMatrixSplice::CollectionMatrixSplice(argtype args) :
  CollectionMatrixSplice(&args) {
  check_all_used(args);
}

void CollectionMatrixSplice::set(const int index, std::shared_ptr<MonteCarlo> mc) {
  if (index >= num()) {
    clones_.resize(index + 1);
  }
  clones_[index] = mc;
}

const MonteCarlo& CollectionMatrixSplice::clone(const int index) const {
  ASSERT(index < num(), "index: " << index << " >= num: " << num());
  return const_cast<MonteCarlo&>(*clones_[index]);
}

MonteCarlo * CollectionMatrixSplice::get_clone(const int index) {
  ASSERT(index < num(), "index: " << index << " >= num: " << num());
  return clones_[index].get();
}

FlatHistogram CollectionMatrixSplice::flat_histogram(const int index) const {
  std::stringstream ss;
  clone(index).criteria().serialize(ss);
  return FlatHistogram(ss);
}

void CollectionMatrixSplice::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2974, ostr);
  feasst_serialize(clones_, ostr);
  feasst_serialize(min_window_size_, ostr);
  feasst_serialize(hours_per_, ostr);
  feasst_serialize(ln_prob_file_, ostr);
  feasst_serialize(ln_prob_file_append_, ostr);
  feasst_serialize(checkpoint_, ostr);
  feasst_serialize_endcap("CollectionMatrixSplice", ostr);
}

CollectionMatrixSplice::CollectionMatrixSplice(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2974, "version: " << version);
  // HWH for unknown reasons, this does not work
  //feasst_deserialize(&clones_, istr);
  int dim1;
  istr >> dim1;
  clones_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      clones_[index] = std::make_shared<MonteCarlo>(istr);
    }
  }
  feasst_deserialize(&min_window_size_, istr);
  feasst_deserialize(&hours_per_, istr);
  feasst_deserialize(&ln_prob_file_, istr);
  feasst_deserialize(&ln_prob_file_append_, istr);
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(checkpoint_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      checkpoint_ = std::make_shared<Checkpoint>(istr);
    }
  }
  feasst_deserialize_endcap("CollectionMatrixSplice", istr);
}

void CollectionMatrixSplice::set(std::shared_ptr<Checkpoint> checkpoint) {
  checkpoint_ = checkpoint;
}

bool CollectionMatrixSplice::are_all_complete() const {
  for (auto clone : clones_) {
    if (!clone->criteria().is_complete()) {
      return false;
    }
  }
  return true;
}

void CollectionMatrixSplice::run(const double hours) {
  #ifdef _OPENMP
  #pragma omp parallel
  {
    const int thread = omp_get_thread_num();
    if (thread < num()) {
      MakeRun({{"for_hours", str(hours)}})->run(clones_[thread].get());
    }
  }
  #else // _OPENMP
    for (auto clone : clones_) {
       MakeRun({{"for_hours", str(hours)}})->run(clone.get());
    }
  #endif // _OPENMP
}

void CollectionMatrixSplice::write(const std::string& file_name) const {
  if (!file_name.empty()) {
    std::ofstream file;
    if (ln_prob_file_append_) {
      file.open(file_name, std::ofstream::out | std::ofstream::app);
      file << "#\"cpu_hours\":" << cpu_hours() << std::endl;
    } else {
      file.open(file_name, std::ofstream::out);
    }
    if (file.good()) {
      auto tm = MakeTransitionMatrix({{"min_sweeps", "0"}});
      tm->set_cm(collection_matrix());
      file << "state," << tm->write_per_bin_header() << std::endl;
      for (int bin = 0; bin < static_cast<int>(tm->collection().matrix().size()); ++bin) {
        file << bin << "," << tm->write_per_bin(bin) << std::endl;
      }
    }
  }
}

void CollectionMatrixSplice::run_until_all_are_complete() {
  #ifdef _OPENMP
  #pragma omp parallel
  {
    const int thread = omp_get_thread_num();
    while (!are_all_complete()) {
      if (thread < num()) {
        MakeRun({{"for_hours", str(hours_per_)}})->run(clones_[thread].get());
      }

      #pragma omp barrier
      if (thread == 0) {
        checkpoint_->check(*this);

        // write ln_prob
        write(ln_prob_file_);
        adjust_bounds();
      }
      #pragma omp barrier
    }
  }
  #else // _OPENMP
    FATAL("OMP required for CollectionMatrixSplice::run_until_all_are_complete()");
  #endif // _OPENMP
}

CollectionMatrix CollectionMatrixSplice::collection_matrix(
    const int index) const {
  FlatHistogram fh = flat_histogram(index);
  try {
    std::stringstream ss;
    fh.bias().serialize(ss);
    TransitionMatrix tm(ss);
    return tm.collection();
  } catch (...) {
    std::stringstream ss;
    fh.bias().serialize(ss);
    WLTM tm(ss);
    return tm.transition_matrix().collection();
  }
}

CollectionMatrix CollectionMatrixSplice::collection_matrix() const {
  CollectionMatrix cm = collection_matrix(0);
  std::vector<std::vector<Accumulator> > data = cm.matrix();
  FlatHistogram fh0 = flat_histogram(0);
  //ASSERT(fh0.macrostate().soft_min() == 0, "soft min should be 0");
  int last_max = fh0.macrostate().soft_max();
  for (int cli = 1; cli < num(); ++cli) {
    FlatHistogram fh = flat_histogram(cli);
    ASSERT(fh.macrostate().soft_min() == last_max + 1,
      "soft_min: " << fh.macrostate().soft_min() <<
      " last_max: " << last_max);
    last_max = fh.macrostate().soft_max();
    CollectionMatrix cmi = collection_matrix(cli);
    int max_bin = fh.macrostate().soft_max();
    if (cli == num() - 1) {
      max_bin = fh.macrostate().histogram().size() - 1;
    }
    for (int bin = fh.macrostate().soft_min(); bin <= max_bin; ++bin) {
      data[bin] = cmi.matrix()[bin];
    }
  }
  return CollectionMatrix(data);
}

LnProbability CollectionMatrixSplice::ln_prob() const {
  CollectionMatrix cm = collection_matrix();
  LnProbability ln_prob;
  ln_prob.resize(cm.matrix().size());
  cm.compute_ln_prob(&ln_prob);
  return ln_prob;
}

// HWH enable adjusting bounds of a single simulation
void CollectionMatrixSplice::adjust_bounds() {
  if (min_window_size_ <= 0) return;
  bool left_most = true;
  bool right_most = false;
  if (num() == 1) {
    clones_[0]->adjust_bounds(true, true, min_window_size_, NULL);
  } else {
    for (int sim = 0; sim < num() - 1; ++sim) {
      DEBUG("sim " << sim);
      if (sim == num() - 2) {
        right_most = true;
      }
      clones_[sim]->adjust_bounds(left_most, right_most, min_window_size_,
                                  clones_[sim + 1].get());
      left_most = false;
    }
  }
}

}  // namespace feasst
