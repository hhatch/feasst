#include <stdio.h>
#include "feasst/include/feasst.h"

using namespace feasst;

int main() {
  std::cout << "version: " << feasst::version() << std::endl;
  std::string line;
  arglist list;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] != '#') {
      std::stringstream ss(line);
      std::string major;
      ss >> major;
      argtype args;
      while(!ss.eof()) {
        std::string minor, value;
        ss >> minor >> value;
        args[minor] = value;
      }
      list.push_back(std::pair<std::string, argtype>(major, args));
    }
  }
  auto mc = feasst::MakeMonteCarlo(list);
}
