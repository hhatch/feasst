#include "utils/test/utils.h"
#include "models/include/square_well.h"

namespace feasst {

TEST(SquareWell, serialize) {
  SquareWell model;
  std::shared_ptr<Model> model2 = test_serialize<SquareWell, Model>(model,
    "SquareWell 2094 -1 -1 -1 -1 553 ");
}

}  // namespace feasst
