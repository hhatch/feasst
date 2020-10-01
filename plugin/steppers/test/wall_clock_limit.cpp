#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/seek_analyze.h"

namespace feasst {

TEST(WallClockLimit, limit) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
  TRY(
    mc.add(MakeWallClockLimit({{"max_hours", "1e-9"}}));
    mc.attempt(1e4);
    CATCH_PHRASE("exceed the maximum");
  );
  auto limit = test_serialize(*MakeWallClockLimit({{"max_hours", "1e-9"}}));
  EXPECT_EQ(0, SeekAnalyze().index("Log", mc)[0]);
  EXPECT_EQ(0, SeekAnalyze().index("Log", mc)[1]);
  EXPECT_EQ(0, SeekAnalyze().index("Movie", mc)[0]);
  EXPECT_EQ(1, SeekAnalyze().index("Movie", mc)[1]);
  EXPECT_EQ(1, SeekAnalyze().index("WallClockLimit", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("WallClockLimit", mc)[1]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", mc)[0]);
}

}  // namespace feasst
