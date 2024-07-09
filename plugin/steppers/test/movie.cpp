#include "utils/test/utils.h"
#include "configuration/include/file_vmd.h"
#include "configuration/include/file_xyz.h"
#include "steppers/include/movie.h"

namespace feasst {

TEST(Movie, serialize) {
  auto movie = MakeMovie({{"output_file", "tmp"}});
  auto movie2 = test_serialize(*movie);
  //auto movie2 = test_serialize<Movie, Analyze>(*movie);
}

}  // namespace feasst
