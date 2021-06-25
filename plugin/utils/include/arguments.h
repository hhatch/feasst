
#ifndef FEASST_UTILS_ARGUMENTS_H_
#define FEASST_UTILS_ARGUMENTS_H_

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "utils/include/io.h"

namespace feasst {

/// Use a map of string pairs as a dictionary for arguments.
typedef std::map<std::string, std::string> argtype;
typedef std::map<std::string, argtype> arglist;

/**
  Classes should take argtype as input for (optional) arguments that may or may
  not have default values:

  MakeTestArgs({{"arg", "value"}});

  All arguments are pairs of strings that may be converted to int, double or
  bool via utils/include/io.h.

  Constructurs should take argtype as object and pointer.
  The point implementation removes arguments as they are processed, while the
  argtype constructor passes the arguments to the other constructor, and then
  checks that all argtype were used.
  Thus, the pointer implementation is used for base classes and containers,
  while the object implementation is for the user.

  Please see utils/test/arguments.cpp TestArgs and TEST(Arguments, args)
  for a working example of how to implement arguments in a class.
 */

/// Return true if key is used in args (templated for arglist or argtype).
template <typename T>
bool used(const std::string& key, const T& args) {
  const auto pair = args.find(key);
  if (pair != args.end()) {
    return true;
  }
  return false;
}

/**
  Read an argument, but do not remove it.
  WARNING: This should typically only be used for error checking.
 */
std::string str(const std::string& key, const argtype& args);

/// Read an argument and remove it
std::string str(const std::string& key, argtype * args);

//// Depreciate
//argtype get(const std::string& key, arglist * args);

/// If args contains derived class of T, return factory pointer and remove from
/// args.
template <class T>
std::shared_ptr<T> parse(T * obj, arglist * args) {
  std::shared_ptr<T> new_obj;
  const auto& map = obj->deserialize_map();
  //for (auto& arg : args) {
  for (const auto& mp : map) {
    auto pair = args->find(mp.first);//map.find(arg.first);
    if (pair != args->end()) {
      new_obj = obj->factory(pair->first, &pair->second);
      args->erase(pair);
      return new_obj;
    }
  }
  return new_obj;
}

/// Same as above, but with a default value should key not be in args.
std::string str(const std::string& key, argtype * args,
  const std::string dflt);

/// Read an argument and remove it from args, then return as double.
double dble(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
double dble(const std::string& key, argtype * args,
  const double dflt);

/// Read an argument and remove it from args, then return as double.
int integer(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
int integer(const std::string& key, argtype * args,
  const int dflt);

/// Read an argument and remove it from args, then return as a boolean.
bool boolean(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
bool boolean(const std::string& key, argtype * args,
  const bool dflt);

/// Append to given key
void append(const std::string& key, argtype * args, const std::string& append);

/// Check that all arguments are used.
void check_all_used(const argtype& args);

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENTS_H_
