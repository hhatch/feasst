
#ifndef FEASST_MATH_FORMULA_POLYNOMIAL_H_
#define FEASST_MATH_FORMULA_POLYNOMIAL_H_

#include <string>
#include <memory>
#include <vector>
#include "math/include/formula.h"

namespace feasst {

/**
  \f$ f(x) = \sum_{i=0}^n A_n(x - x_0)^i \f$
  where A_0 = f(x_0).
 */
class FormulaPolynomial : public Formula {
 public:
  //@{
  /** @name Arguments
    - coeff[i]: add the i-th coefficient where i corresponds with the power of
      the polynomial.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
   */
  explicit FormulaPolynomial(argtype args = argtype());
  explicit FormulaPolynomial(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  FormulaPolynomial& set_A(const int index, const double A) {
    // resize A coefficients if necessary
    if (static_cast<int>(A_.size()) <= index) {
      A_.resize(index + 1);
    }
    A_[index] = A;
    return *this;
  }

  double evaluate(const double x) const override;
  double derivative(const double x) const override;

  std::unique_ptr<Formula> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit FormulaPolynomial(std::istream& istr);
  virtual ~FormulaPolynomial() {}

  //@}
 private:
  std::vector<double> A_;
};

inline std::unique_ptr<FormulaPolynomial> MakeFormulaPolynomial(
    argtype args = argtype()) {
  return std::make_unique<FormulaPolynomial>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_FORMULA_POLYNOMIAL_H_
