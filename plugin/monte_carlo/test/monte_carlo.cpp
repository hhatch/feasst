#include <memory>
#include <cmath>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/dont_visit_model.h"
#include "system/include/ideal_gas.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_volume.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/test/monte_carlo_benchmark.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/profile_trials.h"
#include "steppers/include/volume.h"
#include "opt_lj/include/visit_model_opt_lj.h"
//#include "cluster/include/energy_map_neighbor.h"
//#include "cluster/include/energy_map_all.h"

namespace feasst {

//TEST(MonteCarlo, tutorial) {
//  MonteCarlo mc(
//    MakeRandomMT19937({{"seed", "time"}}),
//    Configuration(feasst::MakeDomain({{"cubic_box_length", "8"}}),
//      {{"particle_type", install_dir() + "/forcefield/lj.fstprt")}}));
//  mc.add(feasst::Potential(feasst::MakeLennardJones()));
//  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
//  mc.add(feasst::MakeMetropolis(
//    {{"beta", "1.2"}, {"chemical_potential", "1."}}));
//  mc.add(feasst::MakeTrialTranslate(
//    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
//  const int steps_per = 1e3;
//  mc.add(feasst::MakeTune({{"steps_per", feasst::str(steps_per)}}));
//  mc.seek_num_particles(50);
//  mc.add(feasst::MakeLog({{"steps_per", feasst::str(steps_per)}}));
//  mc.add(feasst::MakeMovie(
//   {{"steps_per", feasst::str(steps_per)}, {"file_name", "movie.xyz"}}));
//  mc.add(feasst::MakeCheckEnergy(
//   {{"steps_per", feasst::str(steps_per)}, {"tolerance", "1e-8"}}));
//  mc.attempt(1e5);
//}

TEST(MonteCarlo, serialize) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.set(MakeCheckpoint({{"num_hours", "0.0001"}, {"file_name", "tmp/ljrst"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyze(0).class_name(), "AnalyzeFactory");
  EXPECT_EQ(mc2.analyze(0).analyze(0).class_name(), "Log");
  EXPECT_EQ(mc2.analyze(0).analyze(1).class_name(), "Movie");
  EXPECT_EQ(mc2.modify(0).class_name(), "ModifyFactory");
  EXPECT_EQ(mc2.modify(0).modify(0).class_name(), "CheckEnergy");
  EXPECT_EQ(mc2.modify(0).modify(1).class_name(), "Tune");

  auto mc3 = MakeMonteCarlo({{
    {"Configuration", {{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialTransfer", {{"weight", "4."}, {"particle_type", "0"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"file_name", "tmp/ljrst"}}},
    {"Log", {{"steps_per", str(1e4)}, {"file_name", "tmp/lj.txt"}}},
    {"Movie", {{"steps_per", str(1e4)}, {"file_name", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{"steps_per", str(1e4)}}},
  }});

  MonteCarlo mc4 = test_serialize(*mc3);
  EXPECT_EQ(mc4.trial(0).class_name(), "TrialTranslate");
  EXPECT_EQ(mc4.trial(0).weight(), 1.);
  EXPECT_EQ(mc4.trial(1).class_name(), "TrialAdd");
  EXPECT_EQ(mc4.trial(1).weight(), 2.);
  EXPECT_EQ(mc4.trial(2).class_name(), "TrialRemove");
  EXPECT_EQ(mc4.trial(2).weight(), 2.);
  EXPECT_EQ(mc4.analyze(0).class_name(), "Log");
  EXPECT_EQ(mc4.analyze(1).class_name(), "Movie");
  EXPECT_EQ(mc4.modify(0).class_name(), "CheckEnergy");
  EXPECT_EQ(mc4.modify(1).class_name(), "Tune");
}

TEST(MonteCarlo, NVT_NO_FEASST_BENCHMARK_LONG) {
  srand(123);
  RandomBenchmark random;
  serial(&random);
}

TEST(MonteCarlo, NVT_BENCHMARK_LONG) {
  for (const std::string type : {"not"}) {
  //for (const std::string type : {"opt", "not"}) {
  //for (const std::string type : {"neighbor"}) {
  //for (const std::string type : {"all"}) {
  //for (const std::string type : {"opt", "neighbor", "all"}) {
    INFO("type: " << type);
    MonteCarlo mc;
    mc.set(MakeRandomMT19937({{"seed", "123"}}));
    mc.add(MakeConfiguration({{"cubic_box_length", "8"},
      {"particle_type0", "../forcefield/lj.fstprt"},
      {"xyz_file", "../plugin/monte_carlo/test/data/bench.xyz"}}));
    mc.add(MakePotential(MakeLennardJones()));
//    if (type == "all") {
//      mc.set(0, MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
//    } else if (type == "neighbor") {
//      mc.set(0, MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
//    }
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.set(MakeMetropolis());
    INFO(mc.criteria().current_energy());
    mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  //  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  //  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
    //mc.set(0, Potential(MakeLennardJones(),
    //  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    if (type == "opt") {
      mc.add_to_optimized(MakePotential(MakeLennardJones(), //HWH: prevents ModelEmpty... how to remove?
                                        MakeVisitModelOptLJ()));
    }
    mc.attempt(1e6);  // 4.5s with 50 (see opt_lj for 3s)
    //mc.attempt(1e6);  // 5.4s with 50 (see opt_lj for 4.3s)
    // mc.seek_num_particles(450);
    // mc.attempt(1e5);  // 15 sec with 450 on slow computer
    //mc.attempt(1e3);
    // DEBUG("\n" << mc.timer_str());
  }
}

// HWH Removed due to DCCB issue
//TEST(MonteCarlo, NVT_cell_BENCHMARK_LONG) {
//  MonteCarlo mc;
//  mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"}}),
//                       {{"particle_type", "../forcefield/lj.fstprt"}}));
//  mc.add(MakePotential(MakeLennardJones()));
//  mc.add_to_reference(MakePotential(MakeLennardJones(), MakeVisitModelCell({{"min_length", "1"}})));
//  mc.set(MakeRandomMT19937({{"seed", "default"}}));
//  mc.set(MakeThermoParams({{"beta", "1.2"}}));
//  mc.set(MakeMetropolis());
//  mc.add(MakeTrialTranslate({
//    {"weight", "1"},
//    {"reference_index", "0"},
//    {"num_steps", "10"},
//    {"tunable_param", "1"}}));
////  const int steps_per = 1e3;
////  mc.add(MakeCheckEnergy({{"steps_per", str(steps_per)}, {"tolerance", "1e-10"}}));
//  SeeeekNumParticles(50)
//    .with_thermo_params({{"beta", "0.1"}, {"chemical_potential", "10"}})
//    .with_metropolis()
//    .with_trial_add()
//    .run(&mc);
//  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/cell"}}));
//  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
//  mc.attempt(1e5);  // 2 sec with 50 particles, 10 steps on hwhdesk after moved to VisitModelCell
//  //mc.attempt(1e5);  // 3 sec with 50 particles, 10 steps on hwhdesk after site cells integer
//  //mc.attempt(1e5);  // 4.4 sec with 50 particles, 10 steps on hwhdesk after cell group opt
//  //mc.attempt(1e5);  // 4.5 sec with 50 particles, 10 steps on hwhdesk after domain shift opt
//  //mc.attempt(1e5);  // 5.6 sec with 50 particles, 10 steps on hwhdesk
//  //mc.attempt(1e5);  // 5.1 sec with 50 particles, 10 steps on hwhdesk after opt cell_id
//  //mc.attempt(1e3);
//  // mc.attempt(1e6);
//}
//TEST(MonteCarlo, NVT_two_cells_BENCHMARK_LONG) {
TEST(MonteCarlo, NVT_cells_BENCHMARK_LONG) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "12"},
                            {"particle_type", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential({{"Model", "LennardJones"}}));
//  mc.add_to_reference(MakePotential(MakeLennardJones(), MakeVisitModelCell({{"min_length", "1"}})));
  mc.set(MakeThermoParams({{"beta", "0.1"}, {"chemical_potential", "10"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({
    {"weight", "1"},
//    {"reference_index", "0"},
//    {"num_steps", "4"},
    {"tunable_param", "1"}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "200"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(MakeThermoParams({{"beta", "1.2"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/cell"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  mc.add_to_optimized(MakePotential(MakeLennardJones(), MakeVisitModelCell({{"min_length", "3"}})));
  mc.initialize_system();
  mc.attempt(1e5);
}

TEST(MonteCarlo, NVT_SRSW) {
  const int nMol = 500;
  const double rho = 1e-3;
  const double length = std::pow(static_cast<double>(nMol)/rho, 1./3.);
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", str(length)}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", str(nMol)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e3)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e3)}, {"tolerance", str(1e-9)}}));
  Accumulator pe;
  const int num_trials = 1e3;
  for (int trial = 0; trial < num_trials; ++trial) {
    mc.attempt(1);  // ~4 seconds
    pe.accumulate(mc.criteria().current_energy());
  }
  // HWH temperature not set
  DEBUG("pe " << pe.average());
  EXPECT_LE(mc.modify(0).modify(0).accumulator().average(), 1e-13);
}

TEST(MonteCarlo, GCMC) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-3"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"weight", "4."}, {"particle_type", "0"}}));
  EXPECT_NEAR(mc.trial(0).weight(), 1, NEAR_ZERO);
  EXPECT_NEAR(mc.trial(1).weight(), 2, NEAR_ZERO);
  EXPECT_NEAR(mc.trial(2).weight(), 2, NEAR_ZERO);
  mc.add(MakeNumParticles({{"steps_per_write", str(1e5)},
                           {"file_name", "tmp/ljnum.txt"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergy({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  auto profile = MakeProfileTrials({{"steps_per_update", str(1e2)},
    {"steps_per_write", str(1e2)},
    {"append", "true"},
    {"file_name", "tmp/lj_profile.txt"}});
  mc.add(profile);
  //mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  const int trials = 1e4;
  //const int trials = 1e6;
  mc.attempt(trials);
  EXPECT_EQ(mc.trials().num(), 3);
  EXPECT_NEAR(mc.trial(0).num_attempts(), trials/5, trials*0.02);
  EXPECT_NEAR(mc.trial(1).num_attempts(), trials*2/5., trials*0.025);
  EXPECT_NEAR(mc.trial(2).num_attempts(), trials*2/5., trials*0.025);
//  const double sum0 = profile->profile()[0].sum();
//  const double sum1 = profile->profile()[1].sum();
//  const double sum2 = profile->profile()[2].sum();
//  EXPECT_NEAR(2*sum0, sum1, 0.1*sum0);
//  EXPECT_NEAR(sum1, sum2, 0.1*sum0);
}

TEST(MonteCarlo, GCMC_cell) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.run(MakeAddReference({{"cutoff", "1"}, {"use_cell", "true"}}));
  EXPECT_EQ(mc.system().num_references(), 1);
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-6"}}));
  mc.add(MakeTrialTransfer(
    { {"particle_type", "0"},
      {"num_steps", "4"},
      {"reference_index", "0"}}));
  mc.add(MakeNumParticles({{"steps_per_write", str(1e5)},
                           {"file_name", "tmp/ljnum.txt"}}));
  //mc.attempt(1e4);
  for (int i = 0; i < 1e4; ++i) {
    mc.attempt(1);
    mc.system().potential(0).visit_model().check(mc.configuration());
  }
}

TEST(MonteCarlo, ConstrainNumParticles) {
  for (const double minimum : {0, 1}) {
    MonteCarlo mc;
    mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
    mc.add(MakePotential(MakeLennardJones()));
    mc.add(MakePotential(MakeLongRangeCorrections()));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.set(MakeMetropolis());
    mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
    mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
    mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
    mc.get_system()->get_configuration()->add_particle_of_type(0);
    mc.set(MakeThermoParams({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
    mc.set(MakeMetropolis(
      MakeConstrainNumParticles({{"minimum", str(minimum)},
                                 {"maximum", str(minimum+1)}})));
    mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
    const int index = mc.num_analyzers();
    mc.add(MakeNumParticles({{"steps_per_write", "10000"}}));
    mc.attempt(14);
    if (minimum == 0) {
      EXPECT_LE(mc.analyze(index).accumulator().average(), 1);
    } else if (minimum == 1) {
      EXPECT_GE(mc.analyze(index).accumulator().average(), 1);
    }
  }
}

TEST(MonteCarlo, GCMC_binary_tune) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.get_system()->get_configuration()->add_particle_type("../forcefield/lj.fstprt", "2");
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential0", "-6"}, {"chemical_potential1", "-8"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}, {"particle_type", "0"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "2."}, {"particle_type", "1"}}));
  mc.add(MakeTrialTransfer({{"weight", "4."}, {"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"weight", "4."}, {"particle_type", "1"}}));
  const std::string steps_per = str(int(1e2));
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/lj"}}));
  mc.add(MakeTune({{"steps_per", steps_per}}));
  mc.attempt(1e4);
  EXPECT_GT(mc.trial(0).stage(0).perturb().tunable().value(), 2.5);
  EXPECT_GT(mc.trial(1).stage(0).perturb().tunable().value(), 2.5);
}

TEST(MonteCarlo, ideal_gas_pressure_LONG) {
  const int num = 1;
  const double beta = 1.2, pressure = 0.2;
  const double volume = (num + 1)/beta/pressure;
  INFO("expected volume: " << volume);
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "3"},
    {"particle_type", "../forcefield/spce.fstprt"},
    {"add_particles_of_type0", str(num)}}));
  mc.add(MakePotential(MakeDontVisitModel()));
  //mc.add(MakePotential(MakeIdealGas()));
  mc.set(MakeThermoParams({{"beta", str(beta)}, {"pressure", str(pressure)}}));
  mc.set(MakeMetropolis());
//  mc.add(MakeTrialTranslate());
  mc.add(MakeTrialVolume({{"tunable_param", "0.5"}}));
  const std::string steps_per = str(int(1e2));
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/ideal_gas"}}));
  mc.add(MakeTune({{"steps_per", steps_per}}));
  mc.attempt(1e3);
  mc.add(MakeVolume({{"steps_per_write", steps_per},
                     {"file_name", "tmp/ideal_gas_volume"}}));
  mc.attempt(1e6);
  const Accumulator& vol = mc.analyzers()[mc.num_analyzers()-1]->accumulator();
  EXPECT_NEAR(vol.average(), volume, 4*vol.block_stdev());
}

TEST(MonteCarlo, arglist_unrecognized) {
  TRY(
    MakeMonteCarlo({{{"Banana", {{}}}}});
    CATCH_PHRASE("Unrecognized argument: Banana");
  );
  TRY(
    MakeMonteCarlo({{{"Metropolis", {{}}}}});
    CATCH_PHRASE("set System before Criteria");
  );
  TRY(
    MakeMonteCarlo({{{"RandomMT19937", {{"this_is_not", "an_expected_argument"}}}}});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{{"Checkpoint", {{"this_is_not", "an_expected_argument"}}}}});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{{"Configuration", {{"particle_type0", "../forcefield/lj.fstprt"}, {"this_is_not", "an_expected_argument"}}}}});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"particle_type0", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"this_is_not", "an_expected_argument"}}}
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{}}},
      {"Tune", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_box_length", "8"}, {"particle_type", "../forcefield/lj.fstprt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{}}},
      {"Tune", {{}}},
      {"Run", {{"this_is_not", "an_expected_argument"}}},
    }});
    CATCH_PHRASE("unused argument");
  );
}

TEST(MonteCarlo, argslist_order) {
  auto mc = MakeMonteCarlo({{
    {"RandomModulo", {{"seed", "123"}}},
    {"Configuration", {{"cubic_box_length", "8"},
                       {"particle_type0", "../forcefield/lj.fstprt"},
                       {"particle_type1", "../forcefield/atom.fstprt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "0.1"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
  }});
  EXPECT_EQ(1.2, mc->thermo_params().beta());
}

TEST(MonteCarlo, arglist) {
  auto mc = MakeMonteCarlo({{
    {"Checkpoint", {{"file_name", "tmp/lj.fst"}}},
    {"RandomModulo", {{"seed", "123"}}},
    {"Configuration", {{"cubic_box_length", "8"},
                       {"particle_type0", "../forcefield/lj.fstprt"},
                       {"particle_type1", "../forcefield/atom.fstprt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "0.2"},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Log", {{"steps_per", str(1e2)}, {"file_name", "tmp/lj.txt"}}},
    {"Movie", {{"steps_per", str(1e2)}, {"file_name", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"steps_per", str(1e2)}, {"tolerance", "1e-8"}}},
    {"Tune", {{"steps_per", str(1e2)}}},
    {"Run", {{"until_num_particles", "50"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"Run", {{"num_attempts", str(1e3)}}},
    {"RemoveModify", {{"name", "Tune"}}},
    {"Run", {{"num_attempts", str(1e3)}}},
    {"WriteCheckpoint", {{}}},
  }});
  EXPECT_EQ(mc->random().class_name(), "RandomModulo");
  EXPECT_EQ(2, mc->configuration().num_particle_types());
  EXPECT_EQ(2, mc->system().unoptimized().num());
  EXPECT_EQ(1.2, mc->thermo_params().beta());
  EXPECT_EQ(1, mc->trials().num());
  EXPECT_EQ("TrialTranslate", mc->trial(0).class_name());
  //EXPECT_EQ("TrialAdd", mc->trial(1).class_name());
  EXPECT_EQ(1, mc->num_modifiers());
  EXPECT_EQ("CheckEnergy", mc->modify(0).class_name());
  EXPECT_LT(100, mc->trials().num_attempts());
  EXPECT_EQ(50, mc->configuration().num_particles());
  mc->run(MakeRemoveTrial({{"all", "true"}}));
  mc->run(MakeRemoveModify({{"all", "true"}}));
  EXPECT_EQ(0, mc->trials().num());
  EXPECT_EQ(0, mc->num_modifiers());
}

//Note: triclinic wraps don't wrap inside triclinic
//TEST(MonteCarlo, triclinic_random_2d) {
//  MonteCarlo mc;
//  mc.set(MakeRandomMT19937({{"seed", "123"}}));
//  mc.add(MakeConfiguration({
//    {"side_length0", "8"}, {"side_length1", "8"}, {"xy", "4"},
//    {"particle_type0", "../forcefield/atom2d.fstprt"}}));
//  mc.add(MakePotential(MakeIdealGas()));
//  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential0", "10000"}}));
//  mc.set(MakeMetropolis());
////  mc.add(MakeTrialTranslate({{"tunable_param", "1"}, {"weight", "10"}}));
//  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
//  mc.run(MakeRun({{"until_num_particles", str(1e3)}}));
////  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
////  mc.attempt(1e4);
//  for (int part = 0; part < mc.configuration().num_particles(); ++part) {
//    mc.get_system()->get_configuration()->wrap_particle(part);
//  }
//  FileXYZ().write_for_vmd("tmp/random_2d_triclinic.xyz", mc.configuration());
//}

TEST(MonteCarlo, gen_5_spce_in_triclinic) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({
    {"side_length0", "20"}, {"side_length1", "20"}, {"side_length2", "20"},
    {"xy", "4"}, {"yz", "4"}, {"xz", "4"},
    {"particle_type0", "../forcefield/spce.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential0", "-1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "5"}}));
  FileXYZ().write_for_vmd("tmp/spce_triclinic.xyz", mc.configuration());
}

}  // namespace feasst
