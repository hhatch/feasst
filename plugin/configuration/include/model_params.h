
#ifndef FEASST_CONFIGURATION_MODEL_PARAMS_H_
#define FEASST_CONFIGURATION_MODEL_PARAMS_H_

#include <memory>
#include <vector>
#include <string>
#include "configuration/include/physical_constants.h"
#include "configuration/include/particle.h"
#include "configuration/include/properties.h"

namespace feasst {

class ModelParams;

/**
  Model parameters depend upon site types, such as epsilon, sigma, etc.
  These parameters may also be mixed between two different site types.
  Each model parameter has an assumed mixing behavior.
 */
class ModelParam {
 public:
  ModelParam() { set_name("ModelParam"); }

  /// Add a new site type.
  void add(const double value);

  /// Set the value of the site type.
  void set(const int site_type, const double value) {
    values_[site_type] = value; }

  /// Add a new site type. If the site does not contain the model parameter,
  /// then add the default value.
  void add(const Site site, const double default_value = 0);

  /// Add all site types in particle.
  void add(const Particle particle);

  /// Compute the mixed values.
  void mix();

  /// Set the mixed value of the site types.
  void set_mixed(const int site_type1,
    const int site_type2,
    const double value);

  /// Return the value.
  double value(const int type) const;

  /// Return the values.
  const std::vector<double>& values() const { return values_; }

  /// Return the mixed value.
  double mixed_value(const int type1, const int type2) const {
    return mixed_values_[type1][type2];
  }

  /// Return the mixed values.
  const std::vector<std::vector<double> >& mixed_values() const {
    return mixed_values_; }

  /// Return the number of values.
  int size() const { return static_cast<int>(values_.size()); }

  /// Return the name of the model parameter (e.g., epsilon, cutoff)
  std::string name() const { return name_; }

  /// Return the maximum.
  double max() const { return max_value_; }

  /// Return the mixed maximum.
  double mixed_max() const;

  virtual double compute(const int type1, const int type2,
    const ModelParams& model_params) { return 0.; }

  virtual double compute(const int type1, const ModelParams& model_params) {
    return 0.; }

  /// Define new parameters modeled after the existing ones
  virtual void set_param(const ModelParams& existing);

  /// Set the name of the model parameter.
  ModelParam& set_name(const std::string name) { name_ = name; return *this; }

  /// Return as a human readable string.
  std::string str() const;

  void serialize(std::ostream& ostr) const;
  explicit ModelParam(std::istream& istr);
  virtual ~ModelParam() {}

 private:
  std::string name_;
  std::vector<double> values_;
  std::vector<std::vector<double> > mixed_values_;
  double max_value_;
  double max_mixed_value_;
  std::vector<std::vector<bool> > is_mixed_override_;

  /// Define mixing rules in the derived class.
  /// The default is a simple average, unless one of the values is zero.
  virtual double mix_(const double value1, const double value2);

  // resize is_mixed_override_ by setting to false by default but not
  // overwriting any that were previously set to true.
  void override_resize_();
};

/**
 The epsilon parameter is named "epsilon" in LMP-like data file Site Properties.
 The epsilon parameter has the default mixing rule:
 \f$ \epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j} \f$
 */
class Epsilon : public ModelParam {
 public:
  Epsilon() { set_name("epsilon"); }
  explicit Epsilon(std::istream& istr) : ModelParam(istr) {}

 private:
  double mix_(const double value1, const double value2) override;
};

/**
  The sigma parameter is named "sigma" in LMP-like data file Site Properties.
  The sigma parameter has the default mixing rule:

  \f$ \sigma_{ij} = \left\{
    \begin{array}{lr}
      0 & : \sigma_i \sigma_j = 0 \\
      0.5(\sigma_i + \sigma_j) & : \sigma_i\sigma_j \neq 0
    \end{array}
  \right. \f$
 */
class Sigma : public ModelParam {
 public:
  Sigma() { set_name("sigma"); }
  explicit Sigma(std::istream& istr) : ModelParam(istr) {}
};

/**
  The cut off parameter is named "cutoff" in LMP-like data file Site Properties.
  The cut off parameter has the default mixing rule:

  \f$ r^c_{ij} = \left\{
    \begin{array}{lr}
      0 & : r^c_i r^c_j = 0 \\
      0.5(r^c_i + r^c_j) & : r^c_ir^c_j \neq 0
    \end{array}
  \right. \f$
 */
class CutOff : public ModelParam {
 public:
  CutOff() { set_name("cutoff"); }
  explicit CutOff(std::istream& istr) : ModelParam(istr) {}
};

/**
 The charge parameter is named "charge" in LMP-like data file Site Properties.
 The charge parameter, q, has the default mixing rule:
 \f$ q_{ij} = q_i q_j \f$
 */
class Charge : public ModelParam {
 public:
  Charge() { set_name("charge"); }
  explicit Charge(std::istream& istr) : ModelParam(istr) {}

 private:
  double mix_(const double value1, const double value2) override {
    return value1*value2; }
};

/**
  Container for all model parameters.
 */
class ModelParams : public PropertiedEntity {
 public:
  ModelParams();

  /// Deep copy constructor
  ModelParams(const ModelParams& params);

  /// Add all properties in site.
  void add(const Site site);

  /// Add all site types in particle.
  void add(const Particle particle);

  /// Compute the mixed values.
  void mix();

  /// Return the number of values.
  int size() const;

  /// Modify model parameter of a given site type and name to value.
  void set(const std::string name, const int site_type, const double value);

  /// Modify the mixed model parameter of a pair of given site types and name
  /// to value.
  void set(const std::string name, const int site_type1, const int site_type2,
    const double value);

  /// Return model parameters of specific types.
  const Epsilon& epsilon() const { return const_cast<Epsilon&>(*epsilon_); }
  const Sigma& sigma() const { return const_cast<Sigma&>(*sigma_); }
  const CutOff& cutoff() const { return const_cast<CutOff&>(*cutoff_); }
  const Charge& charge() const { return const_cast<Charge&>(*charge_); }

  /// Return constant pointers for optimized model inner loops.
  const std::vector<std::vector<double> >& mixed_epsilon() const {
    return epsilon_->mixed_values(); }
  const std::vector<std::vector<double> >& mixed_sigma() const {
    return sigma_->mixed_values(); }
  const std::vector<std::vector<double> >& mixed_cutoff() const {
    return cutoff_->mixed_values(); }
  const std::vector<std::vector<double> >& mixed_charge() const {
    return charge_->mixed_values(); }

  /// Add a custom model parameter
  void add(std::shared_ptr<ModelParam> param) {
    params_.push_back(param); }

  /// Return the model parameter with the corresponding name.
  const std::shared_ptr<ModelParam> select(const std::string name) const;

  /// Set the minimum cutoff to sigma.
  /// This is used for HardSphere potentials that don't assign cutoff.
  void set_cutoff_min_to_sigma();

  /// Set the physical constants.
  void set_physical_constants(std::shared_ptr<PhysicalConstants> constants =
    std::make_shared<CODATA2018>());

  /// Return the physical constants.
  const PhysicalConstants& physical_constants() const {
    return const_cast<PhysicalConstants&>(*physical_constants_); }

  /// Return the physical constants.
  const PhysicalConstants& constants() const { return physical_constants(); }

  /// Check
  void check() const;

  /// Return as a human readable string.
  std::string str() const;

  void serialize(std::ostream& ostr) const;
  explicit ModelParams(std::istream& istr);

 private:
  /// All types of model parameters listed here.
  /// If adding new types, make sure to add them to params_ below.
  std::shared_ptr<Epsilon> epsilon_;
  std::shared_ptr<Sigma> sigma_;
  std::shared_ptr<CutOff> cutoff_;
  std::shared_ptr<Charge> charge_;

  /// Add all of the above to this list.
  std::vector<std::shared_ptr<ModelParam> > params_;

  std::shared_ptr<PhysicalConstants> physical_constants_;

  /// Add built-in types to params
  void add_();

  std::shared_ptr<ModelParam> select_(const std::string name);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_MODEL_PARAMS_H_
