#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "configuration/include/utils.h"
#include "configuration/include/file_xyz.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Configuration, type_to_file_name) {
  Configuration config({
    {"particle_type0", "../forcefield/data.atom"},
    {"particle_type1", "../forcefield/data.lj"},
    {"particle_type2", "../forcefield/data.spce"},
  });
  TRY(
    auto config2 = config;
    config2.add_particle_type("../forcefield/data.atom");
    CATCH_PHRASE("already provided");
  );
  EXPECT_EQ(3, config.num_particle_types());
  EXPECT_EQ("../forcefield/data.atom", config.type_to_file_name(0));
  EXPECT_EQ("../forcefield/data.lj", config.type_to_file_name(1));
  EXPECT_EQ("../forcefield/data.spce", config.type_to_file_name(2));
}

TEST(Configuration, non_cubic_domain) {
  Configuration config;
  config.set(MakeDomain({
    {"side_length0", "3"},
    {"side_length1", "4"},
    {"side_length2", "5"},
  }));
  EXPECT_EQ(config.dimension(), 3);
  EXPECT_EQ(config.domain().volume(), 3*4*5);
}

TEST(Configuration, coordinates_and_wrapping) {
  Configuration config(MakeDomain({{"cubic_box_length", "5"}}),
    {{"particle_type0", "../forcefield/data.atom"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  Position pos;
  pos.set_to_origin_3D();
  pos.set_coord(0, -583);
  pos.set_coord(1, 83.34);
  pos.set_coord(2, 0.005783);
  Select select(0, config.select_particle(0));
  //select.particle(0, config);
  config.displace_particles(select, pos);
  config.wrap_particle(0);
  EXPECT_EQ(config.num_particles(), 2);
  EXPECT_EQ(config.dimension(), 3);

  Configuration config2 = test_serialize(config);

  { const Position& site_position = config2.select_particle(0).site(0).position();
    EXPECT_NEAR(2.,       site_position.coord(0), 10*NEAR_ZERO);
    EXPECT_NEAR(-1.66,    site_position.coord(1), 10*NEAR_ZERO);
    EXPECT_NEAR(0.005783, site_position.coord(2), 10*NEAR_ZERO);
  }

  // test unwrapped
  config2.init_wrap(false);
  config2.displace_particles(select, pos);
  config2.wrap_particle(0);
  Configuration config3 = test_serialize(config2);
  { const Position& site_position = config3.select_particle(0).site(0).position();
    EXPECT_NEAR(2 - 583,       site_position.coord(0), 10*NEAR_ZERO);
    EXPECT_NEAR(-1.66 + 83.34,    site_position.coord(1), 10*NEAR_ZERO);
    EXPECT_NEAR(0.005783 + 0.005783, site_position.coord(2), 10*NEAR_ZERO);
  }
}

TEST(Configuration, particle_types_lj) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/data.lj"}});
  config->check();
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(1, config->particle_types().num_site_types());
  EXPECT_EQ(1, config->particle_types().num_sites());
  EXPECT_EQ(1, config->particle_type(0).num_sites());
  EXPECT_EQ(0, config->particle_type(0).site(0).type());
  EXPECT_EQ(1, config->unique_types().num_site_types());
  EXPECT_EQ(1, config->unique_types().num_sites());
  EXPECT_EQ(0, config->unique_type(0).site(0).type());
  EXPECT_EQ(3., config->unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(1., config->unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(1., config->unique_type(0).site(0).property("sigma"));
  EXPECT_EQ(1, config->max_sites_in_any_particle());
}

TEST(Configuration, particle_types_spce) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/data.spce"}});
  config->check();
  EXPECT_EQ(2, config->particle_types().num_site_types());
  EXPECT_EQ(3, config->particle_types().num_sites());
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(3, config->particle_type(0).num_sites());
  EXPECT_EQ(0, config->particle_type(0).site(0).type());
  EXPECT_EQ(1, config->particle_type(0).site(1).type());
  EXPECT_EQ(1, config->particle_type(0).site(2).type());
  EXPECT_EQ(2, config->unique_types().num_site_types());
  EXPECT_EQ(2, config->unique_types().num_sites());
  EXPECT_EQ(0, config->unique_type(0).site(0).type());
  EXPECT_EQ(1, config->unique_type(0).site(1).type());
  EXPECT_EQ(10., config->unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(0.650169581, config->unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(3.16555789, config->unique_type(0).site(0).property("sigma"));
  EXPECT_EQ(-0.8476, config->unique_type(0).site(0).property("charge"));
  EXPECT_EQ(10., config->unique_type(0).site(1).property("cutoff"));
  EXPECT_EQ(0., config->unique_type(0).site(1).property("epsilon"));
  EXPECT_EQ(0., config->unique_type(0).site(1).property("sigma"));
  EXPECT_EQ(0.4238, config->unique_type(0).site(1).property("charge"));
  EXPECT_EQ(3, config->max_sites_in_any_particle());
}

TEST(Configuration, bonds_spce) {
  Configuration config = spce_sample1();
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_EQ(1., config.unique_type(0).bond(0).property("length"));
  EXPECT_EQ(0.000001, config.unique_type(0).bond(0).property("delta"));
  EXPECT_EQ(2, config.particle_type(0).bond(0).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(0).site(0));
  EXPECT_EQ(1, config.particle_type(0).bond(0).site(1));
  EXPECT_EQ(2, config.particle_type(0).bond(1).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(1).site(0));
  EXPECT_EQ(2, config.particle_type(0).bond(1).site(1));

  // add an actual particle and see if it has more info than necessary
  config.add_particle_of_type(0);
  EXPECT_EQ(0, config.particle(0).num_bonds());
}

TEST(Configuration, group) {
  auto config = MakeConfiguration(MakeDomain({{"cubic_box_length", "7"}}));
  TRY(
    Configuration config_err(*config);
    config_err.add(MakeGroup({{"add_site_type", "0"}}));
    CATCH_PHRASE("add groups after particle types");
  );
  config->add_particle_type("../forcefield/data.spce");
  config->add_particle_type("../forcefield/data.lj");
  TRY(
    Configuration config_err(*config);
    config_err.add_particle_of_type(0);
    config_err.add_particle_type("../forcefield/data.lj");
    CATCH_PHRASE("types cannot be added after particles");
  );
  config->add(MakeGroup({{"add_site_type", "0"}, {"add_particle_type", "0"}}), "O");
  config->add(MakeGroup({{"add_site_type", "0"}, {"add_particle_type", "1"}}), "H");
  config->add(MakeGroup({{"add_site_type", "2"}, {"add_particle_type", "1"}}), "none");
  EXPECT_TRUE(config->group_select(0).group()->has_property("0"));
  EXPECT_TRUE(config->group_select(1).group()->has_property("O"));
  EXPECT_TRUE(config->group_select(2).group()->has_property("H"));
  EXPECT_TRUE(config->group_select(3).group()->has_property("none"));
  EXPECT_EQ(4, config->group_selects().size());
  for (int part = 0; part < 100; ++part) {
    config->add_particle_of_type(0);
  }
  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz",
                 config.get());
  config->add_particle_of_type(1);
  EXPECT_EQ(2, config->num_particle_types());
  EXPECT_EQ(3, config->num_site_types());
  const Select& sel0 = config->group_select(1);
  EXPECT_EQ(100, sel0.num_particles());
  EXPECT_EQ(100, sel0.num_sites());
  for (int index = 0; index < sel0.num_particles(); ++index) {
    EXPECT_EQ(1, sel0.site_indices(index).size());
  }
  const Select& sel1 = config->group_select(2);
  EXPECT_EQ(0, sel1.num_particles());
  EXPECT_EQ(0, sel1.num_sites());
  const Select& sel2 = config->group_select(3);
  EXPECT_EQ(1, sel2.num_particles());
  EXPECT_EQ(1, sel2.num_sites());

//  TrialSelectParticle tselect;
//  tselect.random_particle(config, 1);
//  EXPECT_EQ(1, tselect.mobile().num_particles());
//  EXPECT_EQ(1, tselect.mobile().num_sites());
//  EXPECT_EQ(0, tselect.mobile().site_index(0, 0));
//  EXPECT_GT(100, tselect.mobile().particle_index(0));
//  tselect.random_particle(config, 2);
//  EXPECT_EQ(0, tselect.mobile().num_particles());
//  EXPECT_EQ(0, tselect.mobile().num_sites());
//  tselect.random_particle(config, 3);
//  EXPECT_EQ(1, tselect.mobile().num_particles());
//  EXPECT_EQ(1, tselect.mobile().num_sites());
//  EXPECT_EQ(100, tselect.mobile().particle_index(0));
//
  EXPECT_EQ(301, config->num_sites());
//  config->remove_particles(select);
//  EXPECT_EQ(300, config->num_sites());
}

TEST(Configuration, cells) {
  Configuration config(MakeDomain({{"cubic_box_length", "7"}}),
    {{"particle_type0", "../forcefield/data.spce"}});
  config.add(MakeGroup({{"add_site_type", "0"}}));
  config.add_particle_of_type(0);
  TRY(
    config.particle(0).site(0).property("cell");
    CATCH_PHRASE("not found");
  );
  auto domain = std::make_shared<Domain>(config.domain());
  domain->init_cells(1);
  domain->init_cells(1.4, 1);
//  config.init_cells(1);
//  config.init_cells(1.4, 1);
  config.set(domain);
  EXPECT_EQ("cell0", config.domain().cells()[0].label());
  EXPECT_EQ(config.domain().cells()[0].num_total(), 7*7*7);
  int cell0 = round(7*7*7/2. - 0.5);
  int cell1 = round(5*5*5/2. - 0.5);
  Site site = config.particle(0).site(0);
  EXPECT_EQ(cell0, round(site.property("cell0")));
  EXPECT_EQ(cell1, round(site.property("cell1")));
  std::vector<int> indices = {0};
  EXPECT_EQ(config.domain().cells(0).particles()[0].num_particles(), 0);
  EXPECT_EQ(config.domain().cells(1).particles()[0].num_particles(), 0);
  EXPECT_EQ(config.domain().cells(0).particles()[cell0].num_particles(), 1);
  EXPECT_EQ(config.domain().cells(1).particles()[cell1].num_particles(), 1);
  double tmp;
  EXPECT_TRUE(config.particle(0).site(1).properties().value("cell0", &tmp));
  EXPECT_FALSE(config.particle(0).site(1).properties().value("cell1", &tmp));
  Position trajectory({-3.49, -3.49, -3.49});

  DEBUG("displacing particles");
  Select select(0, config.select_particle(0));
  //select.particle(0, config);
  config.displace_particles(select, trajectory);
  EXPECT_EQ(0, config.particle(0).site(0).property("cell0"));
  site = config.particle(0).site(0);
  EXPECT_EQ(0, round(site.property("cell0")));
  EXPECT_EQ(0, round(site.property("cell1")));
  DEBUG("cell02 " << config.particle(0).site(2).has_property("cell0"));
  DEBUG("pos " << config.particle(0).site(2).position().coord(0));
  EXPECT_EQ(1, round(config.particle(0).site(1).property("cell0")));
  EXPECT_EQ(6, round(config.particle(0).site(2).property("cell0")));
  EXPECT_FALSE(config.particle(0).site(1).has_property("cell1"));
  EXPECT_FALSE(config.particle(0).site(2).has_property("cell1"));
  EXPECT_EQ(config.domain().cells(0).particles()[cell0].num_particles(), 0);
  EXPECT_EQ(config.domain().cells(1).particles()[cell1].num_particles(), 0);
  DEBUG(site.property("cell0"));
  DEBUG(site.property("cell1"));
  EXPECT_NE(cell1, round(site.property("cell1")));

  // serialize
  Configuration config2 = test_serialize(config);
  EXPECT_EQ(6, round(config2.particle(0).site(2).property("cell0")));
  EXPECT_EQ(config2.domain().cells(1).particles()[cell1].num_particles(), 0);

  config.remove_particle(select);
  config.check();
}

TEST(Configuration, select_particle_by_group) {
  Configuration config = spce_sample1();
  config.add(MakeGroup({{"add_site_type", "0"}}));
  Particle part = config.particle(2, 1);
  EXPECT_EQ(1, part.num_sites());
}

TEST(Configuration, physical_constants) {
  auto config = MakeConfiguration({{"physical_constants", "CODATA2010"}});
  EXPECT_EQ(config->model_params().physical_constants().boltzmann_constant(),
            1.3806488E-23);
  auto config2 = MakeConfiguration();
  EXPECT_EQ(config2->model_params().physical_constants().boltzmann_constant(),
            1.380649E-23);
  TRY(
    auto config3 = MakeConfiguration({{"physical_constants", "bananas"}});
    CATCH_PHRASE("The class name \"bananas\" is not recognized");
  );
}

TEST(Configuration, set_type_type) {
  Configuration config = spce_sample1();
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  EXPECT_EQ(1, config.particle_type(0).site(1).type());
  EXPECT_EQ(1, config.particle(0).site(1).type());
  EXPECT_EQ(1, config.particle(1).site(1).type());
  EXPECT_EQ(2, config.particle(1).num_sites_of_type(1));
  config.set_site_type(0, 1, 0);
  EXPECT_EQ(0, config.particle_type(0).site(1).type());
  EXPECT_EQ(0, config.particle(0).site(1).type());
  EXPECT_EQ(0, config.particle(1).site(1).type());
  EXPECT_EQ(1, config.particle(1).num_sites_of_type(1));
}

}  // namespace feasst
