#include <stdio.h>
#include "feasst/include/feasst.h"

using namespace feasst;

int main() {
  std::cout << "version: " << feasst::version() << std::endl;
  std::string line;
  arglist list;
  while (std::getline(std::cin, line)) {
    if (!line.empty() && line[0] != '#') {
      //std::cout << "line: " << line << std::endl;
      std::stringstream ss(line);
      //std::cout << "ss: " << ss.str() << std::endl;
      std::string major;
      ss >> major;
      //std::cout << "major: " << major << std::endl;
      argtype args;
      while(!ss.eof()) {
        std::string minor, value;
        ss >> minor >> value;
        //std::cout << "minor: " << minor << std::endl;
        //std::cout << "value: " << value << std::endl;
        args[minor] = value;
      }
      list.push_back(std::pair<std::string, argtype>(major, args));
    }
  }
  //INFO(str(list));
  auto mc = feasst::MakeMonteCarlo(list);
}
