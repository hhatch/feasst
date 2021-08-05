#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/flat_histogram.h"
#include "ewald/include/utils.h"
#include "ewald/include/charge_screened.h"
#include "cluster/include/energy_map_all.h"
#include "chain/include/trial_grow.h"
#include "chain/include/check_rigid_bonds.h"

namespace feasst {

double energy_av44(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

// HWH add num steps to spce fh test for DCCB diagnosis
MonteCarlo test_spce_avb_grow_fh(std::shared_ptr<Bias> bias,
    const int num_steps = 1,
    bool test = true,
    const int min = 0,
    const int max = 5,
    const int steps_per = 1e3) {
  const bool avb = true;
  INFO(bias->class_name());
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "123"}}));
  argtype spce_args = {{"physical_constants", "CODATA2010"},
                       {"cubic_box_length", "20"},
                       {"alphaL", "5.6"},
                       {"kmax_squared", "38"},
                       //{"table_size", "0"},
                      };
  int ref = -1;
  if (num_steps > 1) {
    //spce_args.insert({"dual_cut", str(10)});
    spce_args.insert({"dual_cut", str(3.16555789)});
    ref = 0;
  }
  mc.set(spce(spce_args));
  if (avb) {
    auto pot = MakePotential(
      MakeModelTwoBodyFactory({MakeLennardJones(), MakeChargeScreened({{"table_size", "0"}})}),
      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))//,
      //{{"table_size", "1e6"}}
    );
    //mc.set(1, pot);
    mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  }
  const double beta = 1/kelvin2kJpermol(525, mc.configuration()); // mol/kJ
  mc.set(MakeThermoParams({{"beta", str(beta)},
     {"chemical_potential", str(-8.14/beta)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    bias);
  mc.set(criteria);
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
  if (avb) {
    //mc.add(MakeNeighborCriteria({{"maximum_distance", "10"}, {"minimum_distance", "3"}, {"site_type0", "0"}, {"site_type1", "0"}, {"potential_index", "1"}}));
    mc.add(MakeTrialGrow(
      {
        {{"transfer", "true"}, {"particle_type", "0"}, {"weight", "10"}, {"site", "0"}},
        //{{"regrow_avb2", "true"}, {"particle_type", "0"}, {"weight", "0.000010"}, {"site", "0"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "0"}},
        {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
      },
      {{"num_steps", "1"}, {"reference_index", "0"}}
    ));
  } else {
    mc.add(MakeTrialTransfer({
      {"particle_type", "0"},
      {"weight", "4"},
      {"reference_index", str(ref)},
      {"num_steps", str(num_steps)}}));
  }
  SeekNumParticles(min).with_thermo_params({{"beta", "1"}, {"chemical_potential", "1"}}).with_metropolis().run(&mc);
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/spce_fh"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeCheckRigidBonds({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/spce_crit.txt"}}));
  auto energy = MakeEnergy({
    {"file_name", "tmp/spce_fh_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}});
  mc.add(energy);
  MonteCarlo mc2 = test_serialize(mc);
  mc2.run_until_complete();

  if (!test) return mc2;

  EXPECT_LE(mc2.system().configuration().num_particles(), 5);

  // known values of lnpi and energy
  const std::vector<std::vector<double> > lnpi_srsw = {
    {-2.7207, 0.015},
    {-1.8523, 0.015},
    {-1.54708, 0.016},
    {-1.51786, 0.015},
    {-1.6479, 0.015},
    {-1.8786, 0.03}};
  const std::vector<std::vector<double> >  en_srsw = {
    {0, 1e-13},
    {-0.0879115, 1.1293158298007674394e-05},
    {-2.326, 0.12},
    {-6.806, 0.24},
    {-13.499, 0.5},
    {-22.27, 1.0}};

  FlatHistogram fh(mc2.criteria());
  const LnProbability& lnpi = fh.bias().ln_prob();
  for (int macro = 0; macro < lnpi.size(); ++macro) {
    EXPECT_NEAR(lnpi.value(macro), lnpi_srsw[macro][0],
      15*lnpi_srsw[macro][1]);
//      if (bias->class_name() == "TransitionMatrix") {
      const double en_std = std::sqrt(std::pow(en_srsw[macro][1], 2) +
        std::pow(energy->energy().block_stdev(), 2));
      EXPECT_NEAR(energy_av44(macro, mc2), en_srsw[macro][0], 15.*en_std);
//      }
  }

  return mc2;
}

TEST(MonteCarlo, spce_fh2_LONG) {
  test_spce_avb_grow_fh(MakeTransitionMatrix({{"min_sweeps", "25"}}));
}

}  // namespace feasst
