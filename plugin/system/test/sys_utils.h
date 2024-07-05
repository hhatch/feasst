
#ifndef FEASST_TEST_SYSTEM_UTILS_H_
#define FEASST_TEST_SYSTEM_UTILS_H_

#include "utils/include/arguments.h"
#include "configuration/include/group.h"
#include "system/include/system.h"
#include "system/include/potential.h"

namespace feasst {

// HWH depreciate
/*
  Return a system with two Lennard-Jones particles.

  args:
  - cubic_side_length: default 6
 */
System two_particle_system(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_TEST_SYSTEM_UTILS_H_
