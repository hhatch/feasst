#include <cmath>
#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "example/include/analyze_example.h"

namespace feasst {

TEST(AnalyzeExample, serialize) {
  auto movie = MakeAnalyzeExample({{"file_name", "tmp"}});
  auto movie2 = test_serialize<AnalyzeExample, Analyze>(*movie);
}

TEST(AnalyzeExample, ideal_gas_fluid_geometric_center_LONG) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../forcefield/atom.fstprt"},
                       {"cubic_box_length", "8"},
                       {"add_particles_of_type0", "100"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"Metropolis", {{"num_trials_per_iteration", "1e2"}}},
    {"TrialTranslate", {{"tunable_param", "4."}}},
    {"AnalyzeExample", {{"trials_per_write", str(1e3)},
                        {"file_name", "tmp/ig_center.csv"},
                        {"start_after_iteration", "1"}}},
    {"Run", {{"num_trials", "1e5"}}},
  }});
  std::stringstream ss;
  mc->analyze(0).serialize(ss);
  AnalyzeExample analyze_example(ss);
  EXPECT_EQ(analyze_example.class_name(), "AnalyzeExample");
  for (int dim = 0; dim < mc->configuration().dimension(); ++dim) {
    const Accumulator& center = analyze_example.geometric_center(dim);
    EXPECT_TRUE(std::abs(center.average()) < 3*center.block_stdev());
  }
}

}  // namespace feasst
