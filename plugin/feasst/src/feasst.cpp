#include <stdio.h>
#include "feasst/include/feasst.h"

int main() {
  std::cout << "version: " << feasst::version() << std::endl;
  std::string line;
  feasst::arglist list;
  feasst::argtype variables;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] != '#') {
      std::stringstream ss(line);
      std::string major, minor, value;
      ss >> major;
      feasst::argtype args;
      bool assign_to_list = true;
      while(!ss.eof()) {
        ss >> minor >> value;
        DEBUG("major " << major << " minor " << minor << " value " << value);
        if (major == "set_variable") {
          DEBUG("setting variable");
          variables[minor] = value;
          assign_to_list = false;
        } else if (variables.count(value) > 0) {
          DEBUG("using variable");
          args[minor] = variables[value];
        } else {
          DEBUG("no variable");
          args[minor] = value;
        }
      }
      if (assign_to_list) {
        list.push_back(std::pair<std::string, feasst::argtype>(major, args));
      }
    }
  }
  auto mc = std::make_shared<feasst::MonteCarlo>(list);
  return 0;
}
