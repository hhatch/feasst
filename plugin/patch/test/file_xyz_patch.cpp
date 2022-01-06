#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "system/include/hard_sphere.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/log_and_movie.h"
#include "models/include/square_well.h"
#include "mayer/include/mayer_sampling.h"
#include "patch/include/visit_model_inner_patch.h"
#include "patch/include/file_xyz_patch.h"

namespace feasst {

TEST(FileXYZPatch, serialize) {
  auto patch = MakeFileXYZPatch({{"append", "true"}});
  FileXYZPatch patch2 = test_serialize(*patch);
  EXPECT_TRUE(patch2.append());
}

TEST(FileXYZPatch, patch) {
  MonteCarlo mc;
  const std::string fstprt = install_dir() + "/plugin/patch/forcefield/janus.fstprt";
  { auto config = MakeConfiguration({{"cubic_box_length", "8"}});
    config->add_particle_type(fstprt);
    config->add_particle_type(fstprt, "2");
    config->add(MakeGroup({{"site_type0", "0"}, {"site_type1", "2"}}));
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    mc.add(config);
  }
  EXPECT_EQ(2, mc.configuration().num_particles());
  mc.add(MakePotential(MakeSquareWell(),
                       MakeVisitModel(MakeVisitModelInnerPatch()),
                       {{"group_index", "1"}}));
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc.configuration());
}

}  // namespace feasst
