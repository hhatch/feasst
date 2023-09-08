
#ifndef FEASST_MONTE_CARLO_STEPPER_H_
#define FEASST_MONTE_CARLO_STEPPER_H_

#include <string>
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"

namespace feasst {

class Criteria;
class System;
class TrialFactory;

/**
  Perform an action (update or write) every so many trials.
  This action could be read-only (see Analyze) or not (see Modify).
  Write to screen if file name is not provided.
 */
class Stepper {
 public:
  /**
    args:
    - trials_per_write: Set the number of trials per write (default: 1).
      Disabled if negative value is provided.
    - trials_per_update: Set the number of trials per update (default: 1).
      Disabled if negative value is provided.
    - file_name: Set the file name to write output (default: empty).
      If empty, write to screen.
    - append: append file output if set to true.
      Do not append if false (default: "false").
    - clear_file: set true to clear contents of file_name, if exists.
      (default: false).
    - stop_after_phase: stop when simulation reaches this phase index.
      If -1, never stop (default: -1).
    - start_after_phase: start when simulation reaches this phase index.
      If -1, start at beginning (default: -1).
    - file_name_append_phase: append phase to file name (default: false)
    - multistate: set "true" to copy for each state (default: false)
    - multistate_aggregate: aggregate the writing of all states, only when
      multistate is enabled (default: true).
      Thus, trials_per_write refers now to the writing of all states.
      Individual states no longer write.
    - stop_after_iteration: stop when Criteria reaches this iteration.
      If -1, never stop (default: -1).
    - start_after_iteration: start when Criteria reaches this iteration.
      If -1, start at beginning (default: -1).
    - rewrite_header: set true to rewrite header every write (default: true).
      If multistate_aggregate, automatically set to false.
    - Accumulator arguments.
    - configuration_index: index of configuration (default: 0)
   */
  explicit Stepper(argtype args = argtype());
  explicit Stepper(argtype * args);

  /// Return the number of trials per update
  int trials_per_update() const { return trials_per_update_; }

  /// Return the number of trials per write.
  int trials_per_write() const { return trials_per_write_; }

  /// Return the file name.
  const std::string file_name() const { return file_name_; }

  /// Return the file name with optionally appended phase.
  std::string file_name(const Criteria& criteria) const;

  /// Return true if phase is to be appended to file name.
  bool file_name_append_phase() const { return file_name_append_phase_; }

  /// Empty the file name.
  void empty_file_name() { file_name_ = ""; }

  /// Return true if appending.
  bool append() const { return append_; }

  /// Return true if header is to be rewritten.
  bool rewrite_header() const { return rewrite_header_; }

  /// Stop after simulation reaches this phase index.
  int stop_after_phase() const { return stop_after_phase_; }

  /// Stop after simulation reaches this phase index.
  int start_after_phase() const { return start_after_phase_; }

  /// Stop after Criteria reaches this iteration.
  int stop_after_iteration() const { return stop_after_iteration_; }

  /// Stop after Criteria reaches this iteration.
  int start_after_iteration() const { return start_after_iteration_; }

  /// Return the configuration index
  int configuration_index() const { return configuration_index_; }

  /// Given the system, return the configuration.
  const Configuration& configuration(const System& system) const {
    return system.configuration(configuration_index_); }

  /// Set the state. Append file name if not empty.
  void set_state(const int state = 0);

  /// Return if multistate.
  bool is_multistate() const { return is_multistate_; }

  /// Return the state.
  int state() const { return state_; }

  /// Return the accumulator.
  const Accumulator& accumulator() const { return accumulator_; }

  /// Return the number of trials since update.
  int trials_since_update() const { return trials_since_update_; }

  /// Return the number of trials since write.
  int trials_since_write() const { return trials_since_write_; }

  /// Return true if aggregating the write of multistate.
  bool is_multistate_aggregate() const { return is_multistate_aggregate_; }

  /// Return the header for writing.
  virtual std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const { return std::string(""); }

  /// Write to standard output if file name is not set. Otherwise, output file.
  void printer(const std::string output, const std::string& file_name);

  virtual std::string class_name() const { return std::string("Stepper"); }

  void serialize(std::ostream& ostr) const;
  Stepper(std::istream& istr);
  virtual ~Stepper() {}

 protected:
  int trials_since_update_ = 0;
  int trials_since_write_ = 0;
  Accumulator accumulator_;

  /// Note that this should not be called after set_state, which appends name.
  void set_file_name(const std::string file_name) { file_name_ = file_name; }

  virtual void set_trials_per_update(const int trials = 1) {
    trials_per_update_ = trials; }

  virtual void set_trials_per_write(const int trials) {
    trials_per_write_ = trials; }

  /// Check if it is time to update or write. Return true if so.
  bool is_time(const int trials_per, int * trials_since);

  /// Set file output to append.
  void set_append() { append_ = true; }

  /// Set file output to not append.
  void set_no_append() { append_ = false; }

  /// Replicate the stepper individually for each state during initialization
  /// of the factory.
  void set_multistate(const bool multi) { is_multistate_ = multi; }

 private:
  int trials_per_update_;
  int trials_per_write_;
  std::string file_name_;
  bool append_;
  int stop_after_phase_;
  int start_after_phase_;
  int stop_after_iteration_;
  int start_after_iteration_;
  bool file_name_append_phase_;
  bool is_multistate_;
  bool is_multistate_aggregate_;
  int state_ = 0;
  int configuration_index_;
  bool rewrite_header_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_STEPPER_H_
