
#ifndef FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_
#define FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_

#include "utils/include/arguments.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Represent anisotropic sites using a tabular potential.
 */
class VisitModelInnerTable : public VisitModelInner {
 public:
  explicit VisitModelInnerTable(argtype args = argtype());
  explicit VisitModelInnerTable(argtype * args);
  void precompute(Configuration * config) override;
  void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<VisitModelInnerTable>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerTable(std::istream& istr);
  virtual ~VisitModelInnerTable() {}
};

inline std::shared_ptr<VisitModelInnerTable> MakeVisitModelInnerTable(
    argtype args = argtype()) {
  return std::make_shared<VisitModelInnerTable>(args);
}

}  // namespace feasst

#endif  // FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_
