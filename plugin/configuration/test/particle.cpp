#include "utils/test/utils.h"
#include "configuration/include/file_lmp.h"

namespace feasst {

TEST(Particle, getset) {
  Position position;
  Site site;
  std::vector<double> x = {3.5, 796.4, -45.4};
  position.set_vector(x);
  site.set_position(position);
}

TEST(Particle, check) {
  Particle particle;
  particle.add(Site(Position({0, 0, 0})));
  Site site;
  Position pos;
  pos.set_vector({0, 0});
  site.set_position(pos);
  TRY(
    particle.add(site);
    particle.check();
    CATCH_PHRASE("size error");
  );
}

TEST(Particle, bond) {
  Particle chain = FileLMP().read("../forcefield/data.chain10");
  TRY(
    chain.bond(9, 10);
    CATCH_PHRASE("not found");
  );
  EXPECT_EQ(0, chain.bond(9, 8).type());
  EXPECT_EQ(8, chain.bond_neighbors(9)[0]);
  EXPECT_EQ(9, chain.bond_neighbors(8)[1]);
  EXPECT_EQ(7, chain.bond_neighbors(8)[0]);
  EXPECT_EQ(0, chain.bond_neighbors(1)[0]);
  EXPECT_EQ(2, chain.bond_neighbors(1)[1]);

  // serialize
  Particle chain2 = test_serialize(chain);
  EXPECT_EQ(10, chain2.num_sites());
  EXPECT_EQ(9, chain2.num_bonds());
  EXPECT_EQ(0, chain2.bond(9, 8).type());
}

TEST(Particle, angle) {
  Particle spce = FileLMP().read("../forcefield/data.spce");
  EXPECT_EQ(0, spce.angle(1, 0, 2).type());
  EXPECT_EQ(0, spce.angle(1, 2, 0).type());
  EXPECT_EQ(0, spce.angle(0, 1, 2).type());
  Particle spce2 = test_serialize(spce);
}

}  // namespace feasst
