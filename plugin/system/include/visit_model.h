
#ifndef FEASST_SYSTEM_VISIT_MODEL_H_
#define FEASST_SYSTEM_VISIT_MODEL_H_

#include <memory>
#include <string>
#include <map>
#include <sstream>
#include "math/include/position.h"
#include "system/include/synchronize_data.h"
#include "system/include/visit_model_inner.h"

namespace feasst {

class Domain;
class Configuration;
class ModelOneBody;
class ModelTwoBody;
class ModelThreeBody;

// HWH document and lint VisitModel
/**
  See Model for a description of the compute methods. These are mirrored by
  simply switching the calling object and the first argument
  (.e.g, Model.compute(Visitor, ...) vs Visitor.compute(Model, ...)
 */
class VisitModel {
 public:
  explicit VisitModel(std::shared_ptr<VisitModelInner> inner =
      std::make_shared<VisitModelInner>()) {
    set_inner(inner); }

  void set_inner(const std::shared_ptr<VisitModelInner> inner) {
    inner_ = inner; }

  const VisitModelInner& inner() const {
    return const_cast<VisitModelInner&>(*inner_); }

  virtual void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelThreeBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);

  // If model parameters are not given, then obtain them from the configuration.
  void compute(
      ModelOneBody * model,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelOneBody * model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelTwoBody * model,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelTwoBody * model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelThreeBody * model,
      Configuration * config,
      const int group_index = 0);

  /// Return the energy.
  double energy() const { return energy_; }

  /// Set the energy.
  void set_energy(const double energy) { energy_ = energy; }

  // HWH only used for EnergyMapNeighborMW
  void final_en(const Select& select) { inner_->final_en(select); }
  void zero_energy() {
    energy_ = 0.;
    inner_->set_energy(0.);
  }

  /// Increment the energy.
  void increment_energy(const double energy) { energy_ += energy; }

  /// Test if energy of whole system is consistent with sum of energy
  /// of selection by particles.
  void check_energy(
      Model * model,
      Configuration * config,
      const int group_index = 0);

  virtual void revert(const Select& select) { inner_->revert(select); }
  virtual void finalize(const Select& select) { inner_->finalize(select); }

  virtual void precompute(Configuration * config) {
    inner_->precompute(config); }

  virtual void check() const { inner_->check(); }

  // Synchronize with another object of the same type.
  // Typically used with prefetch.
  void synchronize_(const VisitModel& visit, const Select& perturbed);
  const SynchronizeData& data() const { return data_; }

  /// Change the volume.
  virtual void change_volume(const double delta_volume, const int dimension) {}

  // serialization
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    ostr << class_name_ << " ";
    serialize_visit_model_(ostr);
  }
  virtual std::shared_ptr<VisitModel> create(std::istream& istr) const {
    return std::make_shared<VisitModel>(istr); }
  std::map<std::string, std::shared_ptr<VisitModel> >& deserialize_map();
  std::shared_ptr<VisitModel> deserialize(std::istream& istr);
  explicit VisitModel(std::istream& istr);
  virtual ~VisitModel() {}

 protected:
  std::string class_name_ = "VisitModel";
  void serialize_visit_model_(std::ostream& ostr) const;
  VisitModelInner * get_inner_() const { return inner_.get(); }

  // HWH hacky addition for optimization: also, prep inner for reverting,
  // because this is called at beginning of every pair-wise selection compute
  // optimization to avoid repeated construction of Position.
  Position relative_, pbc_, origin_;
  void init_relative_(const Domain& domain, Position * relative,
                      Position * pbc);

  SynchronizeData data_;

 private:
  double energy_ = 0.;
  std::shared_ptr<VisitModelInner> inner_;
};

inline std::shared_ptr<VisitModel> MakeVisitModel() {
  return std::make_shared<VisitModel>();
}

inline std::shared_ptr<VisitModel> MakeVisitModel(
    std::shared_ptr<VisitModelInner> inner) {
  return std::make_shared<VisitModel>(inner);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_H_
