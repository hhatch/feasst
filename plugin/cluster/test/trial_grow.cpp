#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/metropolis.h"
#include "chain/include/trial_grow.h"
#include "cluster/include/energy_map_neighbor_criteria.h"

namespace feasst {

// Seems there is an issue with TrialGrow transfer_avb spce
TEST(TrialGrow, transfer_avb_spce) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "20"}}), {
      {"particle_type", "../forcefield/data.spce"}});
    config.add_particle_of_type(0);
    system.add(config);
  }
  const Configuration& config = system.configuration();
  auto ncrit = MakeNeighborCriteria({{"maximum_distance", "4"}, {"minimum_distance", "3.2"}, {"site_type0", "0"}, {"site_type1", "0"}, {"potential_index", "0"}});
  system.add(MakePotential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighborCriteria(ncrit)))));
  system.add_to_reference(MakePotential(MakeDontVisitModel()));
  system.energy();
  system.finalize();
  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "123"}});
  //ran = MakeRandomMT19937({{"seed", "1591972002"}});
  ran = MakeRandomMT19937({{"seed", "1628624844"}});
  auto metropolis = MakeMetropolis();
  system.set(MakeThermoParams({
    {"beta", "0.1"},
    {"chemical_potential", "1"}}));
  system.add(ncrit);
  const double vol_av = system.neighbor_criteria(0).volume(config.dimension());

  INFO("vol_av: " << vol_av);
  auto grow = MakeTrialGrow(
    {
      {{"transfer_avb", "true"}, {"particle_type", "0"}, {"weight", "4"}, {"site", "0"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "0"}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
    },
    {{"num_steps", "1"}, {"reference_index", "0"}}
  );
  grow->precompute(metropolis.get(), &system);
  double en_old = metropolis->current_energy();
  bool accepted = grow->attempt(metropolis.get(), &system, 0, ran.get());
  EXPECT_EQ(grow->num(), 2);
  config.check();
  INFO(config.num_particles());
  if (config.num_particles() != 2) return;
  ASSERT(accepted, "er");
  EXPECT_EQ(config.particle(0).type(), 0);
  EXPECT_EQ(config.particle(1).type(), 0);
  INFO(config.particle(0).site(0).position().str());
  INFO(config.particle(0).site(1).position().str());
  INFO(config.particle(0).site(2).position().str());
  INFO(config.particle(1).site(0).position().str());
  INFO(config.particle(1).site(1).position().str());
  INFO(config.particle(1).site(2).position().str());
  double delta = grow->trial(0).accept().energy_new() - en_old;
  INFO("deltaU " << delta);
  INFO(grow->trial(0).accept().ln_metropolis_prob());
  EXPECT_NEAR(grow->trial(0).accept().ln_metropolis_prob(),
    std::log(vol_av/2.)
    -system.thermo_params().beta()*delta
    +system.thermo_params().beta_mu(0),
    1e-14);

//  system.energy();
//  system.finalize();
  INFO("en: " << system.stored_energy());
  INFO("en: " << metropolis->current_energy());
  en_old = metropolis->current_energy();

  INFO("**begin remove test**");

  accepted = grow->attempt(metropolis.get(), &system, 1, ran.get());
  //ASSERT(accepted, "er");
  delta = - grow->trial(1).accept().energy_old();
  INFO("deltaU: " << delta);
  EXPECT_NEAR(grow->trial(1).accept().ln_metropolis_prob(),
    -std::log(vol_av/2.)
    -system.thermo_params().beta()*delta
    -system.thermo_params().beta_mu(0),
    1e-14);

  en_old = metropolis->current_energy();
  if (accepted) return;

//  INFO("**add a third**");
//  accepted = grow->attempt(metropolis.get(), &system, 0, ran.get());
//  INFO("old en " << en_old);
//  INFO("new en " << grow->trial(0).accept().energy_new());
//  delta = grow->trial(0).accept().energy_new() - en_old;
//  INFO("deltaU: " << delta);
//  EXPECT_DOUBLE_EQ(grow->trial(0).accept().ln_metropolis_prob(),
//    std::log((2./3.)*vol_av/2.)
//    -system.thermo_params().beta()*delta
//    +system.thermo_params().beta_mu(0));
}

}  // namespace feasst
