
#ifndef FEASST_CONFINEMENT_MODEL_SQUARE_WELL_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_SQUARE_WELL_SHAPE_H_

#include "system/include/model_one_body.h"
#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  Note that the input shape of this model represents the shape of the cavity.

  The sigma parameter is the hard sphere.
  The cutoff is the range of the attraction.
  The epsilon parameter is the well depth.
 */
class ModelSquareWellShape : public ModelOneBody,
                     public ShapedEntity {
 public:
  ModelSquareWellShape() {} // serialization only

  // Constructor
  ModelSquareWellShape(
    std::shared_ptr<Shape> shape,
    const argtype& args = argtype());

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelSquareWellShape>(istr); }
  explicit ModelSquareWellShape(std::istream& istr);
  virtual ~ModelSquareWellShape() {}

 private:
  const std::string class_name_ = "ModelSquareWellShape";
  double alpha_;
  Arguments args_;
};

inline std::shared_ptr<ModelSquareWellShape> MakeModelSquareWellShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelSquareWellShape>(shape);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_SQUARE_WELL_SHAPE_H_
