#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

TEST(TrialRotate, spce) {
  System system;
  {
    auto config = MakeConfiguration(MakeDomain({{"cubic_box_length", "4"}}),
      {{"particle_type", "../forcefield/data.spce"}});
    config->add_particle_of_type(0);
    system.add(*config);
  }

  system.add(MakePotential(MakeLennardJones()));
  system.set(MakeThermoParams({{"beta", "1.0"}}));
  Metropolis criteria;
  auto rotate = MakeTrialRotate({{"tunable_param", "90."}, {"weight", "0.5"}});
  EXPECT_EQ(0.5, rotate->weight());
  FileXYZ file;
  file.write("tmp/before", system.configuration());
  RandomMT19937 random;
  rotate->attempt(&criteria, &system, &random);
  file.write("tmp/after", system.configuration());

  test_serialize(*rotate);
}

}  // namespace feasst
