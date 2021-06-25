
#ifndef FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
#define FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

// HWH: Implement "gentle" WL where the bias is updated infrequently.
/**
  Wang Landau flat histogram bias.
  https://doi.org/10.1103/PhysRevLett.86.2050
  https://doi.org/10.1063/1.1615966
 */
class WangLandau : public Bias {
 public:
  /**
    args:
    - min_flatness : Number of flatness checks required for completion.
    - add_to_ln_probability : The initial amount to add to the natural log of
      the macrostate probability upon visiting that state (default: 1.0).
    - reduce_ln_probability : Reduce the amount to add to the natural log of the
      macrostate probability by multiplcation of this factor upon reaching a
      sufficiently flat histogram (default: 0.5).
    - flatness_threshold : The visited states histogram is determined to be flat
      when the percentage between minimum visisted states and average reaches
      this threshold (default: 0.8).
    - updates_per_flat_check: Updates per check for flatness (default: 10^2).
    - min_visit_per_macro: The minimum number of visits for each macrostate
      required during flatness check (default: 10^3).
   */
  explicit WangLandau(argtype args = argtype());
  explicit WangLandau(argtype * args);
  int min_flatness() const { return min_flatness_; }
  void update_or_revert(
    const int macrostate_old,
    const int macrostate_new,
    const double ln_metropolis_prob,
    const bool is_accepted,
    const bool revert) override;
  void set_num_iterations(const int flatness) override;
  const LnProbability& ln_prob() const override {
    return ln_prob_; }
  void resize(const Histogram& histogram) override;
  std::string write() const override;
  std::string write_per_bin(const int bin) const override;
  std::string write_per_bin_header() const override;
  void set_ln_prob(const LnProbability& ln_prob) override;
  const int num_flatness() const { return num_flatness_; }
  std::shared_ptr<Bias> create(std::istream& istr) const override {
    return std::make_shared<WangLandau>(istr); }
  std::shared_ptr<Bias> create(argtype * args) const override {
    return std::make_shared<WangLandau>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WangLandau(std::istream& istr);
  virtual ~WangLandau() {}

 private:
  LnProbability ln_prob_;
  double add_to_ln_probability_ = 0;
  double reduce_ln_probability_ = 0;
  double flatness_threshold_ = 0;
  int updates_per_flat_check_;
  int updates_since_flat_check_ = 0;
  int min_visit_per_macro_;

  /// Count of the number of times a state has been visited since the last time
  /// this histogram was reset after it was deemed to be sufficiently flat.
  std::vector<int> visited_states_;

  /// Number of times the visited states histogram was found to be flat.
  int num_flatness_ = 0;
  int min_flatness_ = 0;

  void flatness_check_();
  /// Perform update when the visited states histogram is found to be flat.
  void flatness_update_();
};

inline std::shared_ptr<WangLandau> MakeWangLandau(
    argtype args = argtype()) {
  return std::make_shared<WangLandau>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_WANG_LANDAU_H_
