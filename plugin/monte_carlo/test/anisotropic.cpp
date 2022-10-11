#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(MonteCarlo, VisitModelInnerTable) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_box_length", "8"}, {"particle_type0", "../plugin/aniso/forcefield/aniso_tabular.fstprt"},
      {"add_particles_of_type0", "1"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "VisitModelInnerTable"}}},
    //{"Potential", {{"VisitModelInner", "VisitModelInnerTable"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"file_name", "tmp/aniso.fst"}}},
    {"Log", {{"trials_per", str(1e0)}, {"file_name", "tmp/aniso.txt"}}},
    {"Movie", {{"trials_per", str(1e0)}, {"file_name", "tmp/aniso.xyz"}}},
    {"CheckEnergy", {{"trials_per", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }});
  EXPECT_EQ(1, mc->configuration().num_particles());
//  EXPECT_NEAR(-2.060346185437E+00, mc->configuration().particle(0).site(0).position().coord(0), NEAR_ZERO);
//  EXPECT_TRUE(std::abs(mc->configuration().particle(0).site(0).position().coord(0)-1.077169909511E+00)>1e-8);
}

}  // namespace feasst
