
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites are
  not available.
  Select a single dihedral from three given anchor sites and one mobile site.
  The mobile site is directly bonded to the first anchor site, the second anchor
  site is bonded to first anchor, and the third anchor is bonded to the second.

  anchor3("l") - anchor2("k") - anchor("j") - mobile("i")

  Dihedral angles are defined as follows, as shown in
  https://en.wikipedia.org/wiki/Dihedral_angle and the first convention in
  http://trappe.oit.umn.edu/torsion.html.
  The first half-plane is defined by all three anchors (l,k,j).
  The second half-plane is defined by the mobile site and the first two anchors (i,j,k).
  The normal of the first plane, \f$n_1\f$, is given by

  \f$n_1=r_{kl} \times r_{jk}\f$

  where

  \f$r_{kl} = r_k - r_l\f$.

  The normal of the second plane, \f$n_2\f$, is given by

  \f$n_2=r_{jk} \times r_{ij}\f$

  and the dihedral angle, \f$\phi\f$, is given by

  \f$\cos\phi = \frac{n_1 \cdot n_2}{|n_1||n_2|}\f$.
 */
class TrialSelectDihedral : public TrialSelectAngle {
 public:
  /**
    args:
    - anchor_site3 : index of third anchor site.
   */
  explicit TrialSelectDihedral(argtype args = argtype());
  explicit TrialSelectDihedral(argtype * args);

  /// Same as TrialSelectAngle, but also add the third anchor site, and add
  /// dihedral_type as an anchor property.
  void precompute(System * system) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectDihedral(std::istream& istr);
  virtual ~TrialSelectDihedral() {}

 protected:
  void serialize_trial_select_dihedral_(std::ostream& ostr) const;

 private:
  int anchor_site3_;
};

inline std::shared_ptr<TrialSelectDihedral> MakeTrialSelectDihedral(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectDihedral>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_DIHEDRAL_H_
