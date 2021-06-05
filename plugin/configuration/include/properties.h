
#ifndef FEASST_CONFIGURATION_PROPERTIES_H_
#define FEASST_CONFIGURATION_PROPERTIES_H_

#include <vector>
#include <string>

namespace feasst {

/**
  Manage custom properties (sites, bonds, etc), typically for use in plugins.
  There is a performance cost associated with accessing the properties by name.
  But accessing the values by index should be done carefully.
 */
class Properties {
 public:
  Properties() {}

  /// Add a new value/name combination.
  void add(const std::string name, const double value);

  /// Add a new value and name, or set its value if name exists.
  void add_or_set(const std::string name, const double value);

  /// Set the value of property name.
  void set(const std::string name, const double value);

  /// Return property value by name.
  double value(const std::string name) const;

  /// Return true if property name exists, and its value.
  bool value(const std::string name, double * value) const;

  /// Return true if property name exists, and its value and index.
  bool value(const std::string name, double * value, int * index) const;

  /// Return true if property name exists.
  bool has(const std::string name) const {
    double val;
    return value(name, &val);
  }

  /// Check that the property values and names are consistent.
  void check() const;

  /// Return all property names.
  std::vector<std::string> names() const { return names_; }

  /// Return all property values.
  const std::vector<double>& values() const { return values_; }

  /// Set value of property by index.
  void set_value(const int index, const double value) {
    values_[index] = value; }

  /// Return the properties as a human readable string.
  std::string str() const;

  /// Return the number of properties.
  int num() const { return static_cast<int>(values_.size()); }

  /// Return true if the properties are equivalent within tolerance.
  bool is_equal(const Properties& properties,
                const double tolerance) const;

  /// Same as above, but with a default tolerance of NEAR_ZERO.
  bool is_equal(const Properties& properties) const;

  void serialize(std::ostream& ostr) const;
  explicit Properties(std::istream& istr);
  ~Properties() {}

 private:
  std::vector<double> values_;
  std::vector<std::string> names_;
};

class PropertiedEntity {
 public:
  PropertiedEntity() {}

  /// Add a new property.
  virtual void add_property(const std::string name, const double value) {
    properties_.add(name, value); }

  /// Add a property, or set its value if name already exists.
  void add_or_set_property(const std::string name, const double value) {
    properties_.add_or_set(name, value); }

  /// Set the value of an existing property.
  void set_property(const std::string name, const double value) {
    properties_.set(name, value); }

  /// Return properties.
  const Properties& properties() const { return properties_; }

  /// Return the property value by name.
  double property(const std::string name) const {
    return properties_.value(name); }

  /// Return true if entity has property of name.
  bool has_property(const std::string name) const {
    return properties_.has(name); }

  /// Set value of property by index.
  void set_property(const int index, const double value) {
    properties_.set_value(index, value); }

  /// Set the properties.
  void set_properties(const Properties& properties) {
    properties_ = properties; }

  void serialize(std::ostream& ostr) const { properties_.serialize(ostr); }
  explicit PropertiedEntity(std::istream& istr) {
    properties_ = Properties(istr); }
  virtual ~PropertiedEntity() {}

 private:
  Properties properties_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_PROPERTIES_H_
