
#ifndef FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_

#include "system/include/model_one_body.h"
#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  Note that the input shape of this model represents the shape of the cavity.

  \f$ U(r) = \epsilon \left( \frac{r}{\sigma} \right)^\alpha \f$
 */
class ModelLJShape : public ModelOneBody,
                     public ShapedEntity {
 public:
  // for serialization only
  ModelLJShape() { class_name_ = "ModelLJShape"; }

  // Constructor
  ModelLJShape(
    std::shared_ptr<Shape> shape,
    /**
      alpha: set the exponent (default: 3).
     */
    argtype args = argtype());

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJShape>(istr); }
  explicit ModelLJShape(std::istream& istr);
  virtual ~ModelLJShape() {}

 private:
  double alpha_;
};

inline std::shared_ptr<ModelLJShape> MakeModelLJShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelLJShape>(shape);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
