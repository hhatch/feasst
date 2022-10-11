
#ifndef FEASST_PATCH_MOVIE_PATCH_H_
#define FEASST_PATCH_MOVIE_PATCH_H_

#include "monte_carlo/include/analyze.h"
#include "aniso/include/file_xyz_euler.h"

namespace feasst {

/**
  Write a trajectory of the site positions and Euler angles.
  Does not overwrite existing file by default.
 */
class TrajEuler : public AnalyzeWriteOnly {
 public:
  explicit TrajEuler(argtype args = argtype());
  explicit TrajEuler(argtype * args);

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  /// Write the configuration.
  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("TrajEuler"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<TrajEuler>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<TrajEuler>(args); }
  TrajEuler(std::istream& istr);

 private:
  FileXYZEuler xyz_;
};

inline std::shared_ptr<TrajEuler> MakeTrajEuler(argtype args = argtype()) {
  return std::make_shared<TrajEuler>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_MOVIE_PATCH_H_
