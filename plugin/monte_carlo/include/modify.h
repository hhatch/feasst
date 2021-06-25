
#ifndef FEASST_MONTE_CARLO_MODIFY_H_
#define FEASST_MONTE_CARLO_MODIFY_H_

#include <memory>
#include <string>
#include <map>
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Perform an action every so many steps that may change the system, criteria
  or trials.
 */
class Modify : public Stepper {
 public:
  Modify() : Stepper() {}
  explicit Modify(argtype * args) : Stepper(args) {}

  /// Initialize and precompute before trials.
  virtual void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) {}

  /// Check every trial if action is to be performed.
  virtual void trial(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory);

  /// Perform update action.
  virtual void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory);

  /// Perform write action.
  virtual std::string write(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory);

  // Access to factory of Modify objects.
  virtual const std::vector<std::shared_ptr<Modify> >& modifiers() const;
  virtual const Modify& modify(const int index) const;

  // serialization
  std::string class_name() const override { return std::string("Modify"); }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Modify> create(std::istream& istr) const;
  virtual std::shared_ptr<Modify> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Modify> >& deserialize_map();
  std::shared_ptr<Modify> deserialize(std::istream& istr);
  std::shared_ptr<Modify> factory(const std::string name, argtype * args);
  explicit Modify(std::istream& istr) : Stepper(istr) {}
  virtual ~Modify() {}

  // HWH only used by ModifyFactory
  void check_update_(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory);
};

/**
  This Modify does not perform writes.
 */
class ModifyUpdateOnly : public Modify {
 public:
  /**
    args:
    - steps_per: update every this many steps
   */
  explicit ModifyUpdateOnly(argtype * args);

  void set_steps_per_write(const int steps) override;

  void set_steps_per(const int steps) { set_steps_per_update(steps); }

  explicit ModifyUpdateOnly(std::istream& istr) : Modify(istr) {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MODIFY_H_
