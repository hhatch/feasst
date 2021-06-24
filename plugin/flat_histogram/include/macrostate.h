
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_H_

#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "math/include/histogram.h"
#include "monte_carlo/include/acceptance.h"

namespace feasst {

/**
  The macrostate is the statistical mechanical variable to which we apply bias
  using flat-histogram methods.
  Thus, for example, in the grand canonical ensemble with a single component,
  the macrostate is the number of particles.
  To apply the flat histogram methods, the macrostate is be broken into a
  contiguous series of "bins".
 */
class Macrostate {
 public:
  // HWH consider depreciating this interface?
  /**
    args:
    - soft_max : optionally, set a soft maximum (default: last histogram bin).
      These soft limits may be changed during a simulation.
      Note that this max is an integer bin number.
    - soft_min : minimum as described above (default: same as histogram).
   */
  Macrostate(const Histogram& histogram, argtype args = argtype());
  Macrostate(const Histogram& histogram, argtype * args);

  /**
    A flattened version of the above constructor that takes Histogram arguments.
   */
  explicit Macrostate(argtype args = argtype());

  /**
    Set the bins of the macrostate by providing a Histogram.
    This is required before the macrostate can be used for flat histogram
    methods.
    The histogram only serves to determine the bins, and should not be
    expanded or have values added during the course of the simulation.
   */
  void set(const Histogram histogram) { histogram_ = histogram; }

  /// Return the histogram.
  const Histogram& histogram() const { return histogram_; }

  /// Return the soft maximum.
  const int soft_max() const { return soft_max_; }

  /// Return the soft minimum.
  const int soft_min() const { return soft_min_; }

  // HWH change this to const references
  /// Return the current value of the macrostate.
  virtual double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const = 0;

  /// Return the current bin of the macrostate.
  int bin(const System& system,
      const Criteria& criteria,
      const Acceptance& acceptance) const {
    return histogram_.bin(value(system, criteria, acceptance)); }

  /// Return the value of the bin.
  double value(const int bin) const { return histogram_.center_of_bin(bin); }

  /// Return whether the current system macrostate is within permissible range
  /// given by the input histogram.
  bool is_allowed(const System& system,
                  const Criteria& criteria,
                  const Acceptance& acceptance) const;

  /// Swap the soft bounds with another macrostate.
  void swap_soft_bounds(Macrostate * macrostate);

  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Macrostate> create(std::istream& istr) const;
  virtual std::shared_ptr<Macrostate> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Macrostate> >& deserialize_map();
  std::shared_ptr<Macrostate> deserialize(std::istream& istr);
  std::shared_ptr<Macrostate> factory(const std::string name, argtype * args);
  explicit Macrostate(std::istream& istr);
  virtual ~Macrostate() {}

 protected:
  void serialize_macrostate_(std::ostream& ostr) const;

 private:
  Histogram histogram_;
  int soft_max_;
  int soft_min_;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_H_
