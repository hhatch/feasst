#include "utils/test/utils.h"
#include "confinement/include/zero_potential.h"

namespace feasst {

TEST(ZeroPotential, serialize) {
  auto obj = MakeZeroPotential();
  ZeroPotential obj2 = test_serialize(*obj);
}

}  // namespace feasst
