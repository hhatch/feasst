
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  This is an example template for implementing your own custom Model.
 */
class ModelExample : public ModelTwoBody {
 public:
  ModelExample();

  /// The energy between two site types depends upon the distance between
  /// them and the model parameters.
  /// See existing Models, such as ModelLJ or ModelLJCutShift for inspiration.
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override {
    return 0;
  }

  /// Serialization is the flatening of an object into a stream of characters
  /// which may be saved to file and later deserialized back into the object.
  void serialize(std::ostream& ostr) const override;

  /// Deserialization is the reverse process of the above.
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelExample>(istr); }
  explicit ModelExample(std::istream& istr);
  virtual ~ModelExample() {}
};

inline std::shared_ptr<ModelExample> MakeModelExample() {
  return std::make_shared<ModelExample>();
}

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
