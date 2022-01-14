#include <stdio.h>
#include "feasst/include/feasst.h"

/**
  Usage: ./rst checkpoint_file.txt
 */
int main(int argc, char ** argv) {
  std::cout
    << "*********************************************************" << std::endl
    << "feasst version: " << feasst::version() << std::endl
    << "*********************************************************" << std::endl;
  ASSERT(argc == 2, "unrecognized number of arguments");
  feasst::MonteCarlo mc;
  feasst::MakeCheckpoint({{"file_name", std::string(argv[1])}})->read(&mc);
  mc.resume();
  return 0;
}
