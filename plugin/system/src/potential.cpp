#include <memory>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/table.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/potential.h"
#include "system/include/model_empty.h"
#include "system/include/model_two_body.h"
#include "system/include/model_two_body_table.h"

namespace feasst {

Potential::Potential(argtype * args) {
  group_index_ = integer("group_index", args, 0);
  if (used("cell_index", *args)) {
    ASSERT(group_index_ == 0, "cell_index overrides group_index");
    group_index_ = integer("cell_index", args);
  }
  prevent_cache_ = boolean("prevent_cache", args, false);
  table_size_ = integer("table_size", args, 0);
}

Potential::Potential(std::shared_ptr<Model> model,
                     argtype args) : Potential(&args) {
  model_ = model;
  visit_model_ = std::make_shared<VisitModel>();
  check_all_used(args);
}

Potential::Potential(std::shared_ptr<VisitModel> visit_model,
                     argtype args) : Potential(&args) {
  model_ = std::make_shared<ModelEmpty>();
  visit_model_ = visit_model;
  check_all_used(args);
}

Potential::Potential(
    std::shared_ptr<Model> model,
    std::shared_ptr<VisitModel> visit_model,
    argtype args) : Potential(&args) {
  model_ = model;
  visit_model_ = visit_model;
  check_all_used(args);
}

Potential::Potential(argtype args) : Potential(&args) {
  model_ = ModelTwoBody().factory(str("Model", &args, "ModelEmpty"), &args);
  visit_model_ =
    VisitModel().factory(str("VisitModel", &args, "VisitModel"), &args);
  check_all_used(args);
}

Potential::Potential(std::shared_ptr<BondVisitor> bond_visitor,
    argtype args) : Potential(args) {
  bond_visitor_ = bond_visitor;
}

void Potential::set(const ModelParams& model_params) {
  model_params_override_ = true;
  model_params_ = ModelParams(model_params);
}

void Potential::set_model_param(const char* name,
    const int site_type,
    const double value) {
  ASSERT(model_params_override_, "you must first initialize model params "
    << "before setting them.");
  model_params_.set(name, site_type, value);
}

const ModelParams& Potential::model_params() const {
  ASSERT(model_params_override_, "When model parameters are not overridden, "
    << "you must also provide the configuration as an argument.");
  return model_params_;
}

const ModelParams& Potential::model_params(const Configuration& config) const {
  if (model_params_override_) {
    return model_params_;
  }
  return config.model_params();
}

double Potential::energy(Configuration * config) {
  ASSERT(bond_visitor_ || visit_model_, "visitor must be set.");
  if (prevent_cache_ || !cache_.is_unloading(&stored_energy_)) {
    if (bond_visitor_) {
      bond_visitor_->compute_all(*config, group_index_);
      stored_energy_ = bond_visitor_->energy();
    } else if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, group_index_, config,
                                       visit_model_.get());
    } else {
      stored_energy_ = model_->compute(group_index_, config, visit_model_.get());
    }
    cache_.load(stored_energy_);
  }
  return stored_energy_;
}

double Potential::select_energy(const Select& select, Configuration * config) {
  ASSERT(visit_model_, "visitor must be set.");
  if (prevent_cache_ || !cache_.is_unloading(&stored_energy_)) {
    if (bond_visitor_) {
      bond_visitor_->compute_all(select, *config);
      stored_energy_ = bond_visitor_->energy();
    } else if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, select, group_index_,
                                       config, visit_model_.get());
    } else {
      stored_energy_ = model_->compute(select, group_index_, config,
                                       visit_model_.get());
    }
    cache_.load(stored_energy_);
  }
  INFO("potential: " << stored_energy_);
  return stored_energy_;
}

int Potential::cell_index() const {
  ASSERT(visit_model_->class_name() == "VisitModelCell", "error");
  return group_index();
}

void Potential::precompute(Configuration * config) {
  if (bond_visitor_) return;
  visit_model_->precompute(config);
  const ModelParams& params = model_params(*config);
  model_->precompute(params);
  const double max_cutoff = maximum(params.cutoff().values());
  const double half_min_side = 0.5*config->domain().min_side_length();
  if (max_cutoff - NEAR_ZERO > half_min_side) {
    WARN("The maximum cutoff:" << max_cutoff << " is greater than half the " <<
         "minimum side length: " << half_min_side);
  }

  if (table_size_ > 0) {
    ASSERT(model_->num_body() == 2, "tables are only implemented for two "
      << "body simulations");
    auto table = MakeModelTwoBodyTable();
    table->set(model_params(*config),
      table_size_,
      config->num_site_types(),
      model_);
    model_ = table;
    table_size_ = 0;
  }
}

void Potential::check(const Configuration& config) const {
  visit_model_->check(config);
}

void Potential::serialize(std::ostream& ostr) const {
  feasst_serialize_version(432, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstdr(visit_model_, ostr);
  feasst_serialize_fstdr(model_, ostr);
  feasst_serialize(stored_energy_, ostr);
  feasst_serialize(model_params_override_, ostr);
  if (model_params_override_) {
    feasst_serialize_fstobj(model_params_, ostr);
  }
  feasst_serialize_fstobj(cache_, ostr);
  feasst_serialize(prevent_cache_, ostr);
  feasst_serialize(table_size_, ostr);
}

Potential::Potential(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(432 == version, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  // feasst_deserialize_fstdr(visit_model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      visit_model_ = visit_model_->deserialize(istr);
    }
  }
  // feasst_deserialize_fstdr(model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      model_ = model_->deserialize(istr);
    }
  }
  feasst_deserialize(&stored_energy_, istr);
  feasst_deserialize(&model_params_override_, istr);
  if (model_params_override_) {
    feasst_deserialize_fstobj(&model_params_, istr);
  }
  feasst_deserialize_fstobj(&cache_, istr);
  feasst_deserialize(&prevent_cache_, istr);
  feasst_deserialize(&table_size_, istr);
}

void Potential::set_model_params(const Configuration& config) {
  set(config.model_params());
}

void Potential::synchronize_(const Potential& potential,
    const Select& perturbed) {
  visit_model_->synchronize_(potential.visit_model(), perturbed);
  stored_energy_ = potential.stored_energy_;
}

}  // namespace feasst
