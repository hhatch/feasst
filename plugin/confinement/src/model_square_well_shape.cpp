#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "confinement/include/model_square_well_shape.h"

namespace feasst {

class MapModelSquareWellShape {
 public:
  MapModelSquareWellShape() {
    ModelSquareWellShape().deserialize_map()["ModelSquareWellShape"] =
      std::make_shared<ModelSquareWellShape>();
  }
};

static MapModelSquareWellShape map_model_hard_shape_ = MapModelSquareWellShape();

ModelSquareWellShape::ModelSquareWellShape(std::shared_ptr<Shape> shape,
  const argtype& args) : ModelOneBody(), ShapedEntity(shape) {
  args_.init(args);
  alpha_ = args_.key("alpha").dflt("3").dble();
}

void ModelSquareWellShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(6482, ostr);
  feasst_serialize(alpha_, ostr);
}

ModelSquareWellShape::ModelSquareWellShape(std::istream& istr)
  : ModelOneBody(), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  feasst_deserialize(&alpha_, istr);
  ASSERT(version == 6482, "unrecognized verison: " << version);
}

double ModelSquareWellShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const {
  const int type = site.type();
  const double sigma = model_params.sigma().value(type);
  const double distance = shape()->nearest_distance(wrapped_site);
  if (distance <= sigma) {
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.epsilon().value(type);
  return -epsilon;
}

}  // namespace feasst
