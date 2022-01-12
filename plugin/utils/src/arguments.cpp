#include <sstream>
#include <cmath>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/io.h"

namespace feasst {

std::string str(const std::string& key, const argtype& args) {
  auto pair = args.find(key);
  if (pair != args.end()) {
    return pair->second;
  }
  return std::string();
}

std::string str(const std::string& key, argtype * args) {
  auto pair = args->find(key);
  ASSERT(pair != args->end(), "key(" << key << ") is required for args but " <<
    "not found in " << str(*args));
  const std::string second = pair->second;
  args->erase(pair);
  return second;
}

//argtype get(const std::string& key, arglist * args) {
//  auto pair = args->find(key);
//  ASSERT(pair != args->end(), "key(" << key << ") is required for args but " <<
//    "not found");// << str(*args));
//  const argtype second = pair->second;
//  args->erase(pair);
//  return second;
//}

std::string str(const std::string& key, argtype * args,
    const std::string dflt) {
  std::string return_str;
  auto pair = args->find(key);
  if (pair != args->end()) {
    return_str = pair->second;
    args->erase(pair);
  } else {
    return_str = dflt;
  }
  return return_str;
}

double dble(const std::string& key, argtype * args) {
  return str_to_double(str(key, args));
}

double dble(const std::string& key, argtype * args,
    const double dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_double(return_str);
  } else {
    return dflt;
  }
}

int integer(const std::string& key, argtype * args) {
  return str_to_int(str(key, args));
}

int integer(const std::string& key, argtype * args,
    const int dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_int(return_str);
  } else {
    return dflt;
  }
}

bool boolean(const std::string& key, argtype * args) {
  return str_to_bool(str(key, args));
}

bool boolean(const std::string& key, argtype * args,
    const bool dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_bool(return_str);
  } else {
    return dflt;
  }
}

void append(const std::string& key, argtype * args, const std::string& append) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string second = pair->second;
    args->erase(pair);
    args->insert({key, second + append});
  }
}

void check_all_used(const argtype& args) {
  if (args.size() != 0) {
    DEBUG(args.size());
    DEBUG(args.begin()->first);
    DEBUG(args.begin()->second);
    ASSERT(args.size() == 1 && args.begin()->first.empty() &&
      args.begin()->second.empty(),
      "unused argument(s): " << str(args) << ". If the arguments are unused " <<
      "then that means the objects did not expect the first keyword " <<
      "supplied in the above argument pair(s). Thus, there was likely a " <<
      "typo or the keyword is intended for a different class.");
  }
}

std::string str(const arglist& args) {
  std::stringstream out;
  out << "{{";
  for (const auto& arg : args) {
    out << "{\"" << arg.first << "\",";
    out << str(arg.second) << "},";
  }
  out << "}}";
  return out.str();
}

void add_if_not_used(const std::string& key, argtype * args,
  const std::string& value) {
  if (!used(key, *args)) {
    args->insert({key, value});
  }
}

argtype extract_until(const std::string key, argtype *args) {
  argtype extracted;
  int num_args = static_cast<int>(args->size());
  for (int iarg = 0; iarg < num_args; ++iarg) {
    auto pair = args->begin();
    if (pair->first == key) {
      args->erase(pair);
      return extracted;
    } else {
      extracted.insert(*pair);
      args->erase(pair);
    }
  }
  return extracted;
}

}  // namespace feasst
