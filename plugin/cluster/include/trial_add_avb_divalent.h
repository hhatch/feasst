
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
Attempt to add a particle of type "t" with site_index_t anywhere in the Domain.
Then add a second particle of type "a" with site_index_a in the AV of
site_index_t.
Finally, add a third particle of type "b" with site_index_a in the AV of
site_index_t.

The derivation of the acceptance criteria follows a similar procedure as
described in TrialAddAVB.

\rst
+-------------------------------------+----------------------------------------+
|Forward event                        |Probability, :math:`\pi_{on}`           |
|                                     |                                        |
|[reverse event]                      |[reverse probability, :math:`\pi_{no}`] |
+-------------------------------------+----------------------------------------+
|Select insert trial                  |:math:`1/w`                             |
|                                     |                                        |
|[Select remove trial]                |:math:`[1/w]`                           |
+-------------------------------------+----------------------------------------+
|Insert site_index_t in box           |:math:`1/V`                             |
|                                     |                                        |
|[Remove site_index_t]                |:math:`1/(N_t + 1 + \delta_{ta} +       |
|                                     |\delta_{tb})`                           |
+-------------------------------------+----------------------------------------+
|Insert site_index_a in AV of         |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Remove site_index_a in AV of        |:math:`1/(N^{s,AV}_a + 1 + \delta_{ab})`|
|site_index_t                         |                                        |
+-------------------------------------+----------------------------------------+
|Insert site_index_b in AV of         |:math:`1/v_{AV}`                        |
|site_index_t                         |                                        |
|                                     |                                        |
|[Remove site_index_b in the AV of    |:math:`1/(N^{s,AV}_b + 1)`              |
|site_index_t]                        |                                        |
+-------------------------------------+----------------------------------------+
|Accept                               |:math:`min(1, \chi)`                    |
|                                     |                                        |
|[Accept]                             |:math:`[min(1, 1/\chi)]`                |
+-------------------------------------+----------------------------------------+

where :math:`N_t` is the number of particles of type t,
:math:`v_{AV}` is the aggregation volume,
:math:`N^{s,AV}_a` is the number of sites with site_index_a, in particles of
type a, that are in the AV of site_index_t and
:math:`\delta_{ta}` is the Kronecker delta as a function of the type index
of the target and added particles (i.e., if add and target are of same type,
then :math:`\delta_{ta} = 1`, otherwise :math:`\delta_{ta} = 0`).

Application of local detailed balance yields the acceptance probability

:math:`\chi =
\frac{V}{N_t + 1 + \delta_{ta} + \delta_{tb}}
\frac{v_{AV}}{N^{s,AV}_a + 1 + \delta_{ab}}
\frac{v_{AV}}{N^{s,AV}_a + 1}
\frac{1}{\Lambda^{3d}}
e^{-\beta\Delta U + \beta\mu_a + \beta\mu_b + \beta\mu_t}`

\endrst
 */
class TrialAddAVBDivalent : public Trial {
 public:
  //@{
  /** @name Arguments
    - particle_type_a: type of second added particle in AV of first.
    - site_index_a: index of site in type a that defines AV (default: 0).
    - particle_type_b: type of third added particle in AV of first.
    - site_index_b: index of site in type b that defines AV (default: 0).
    - SelectParticleAVBDivalent arguments.
   */
  explicit TrialAddAVBDivalent(argtype args = argtype());
  explicit TrialAddAVBDivalent(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddAVBDivalent>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddAVBDivalent>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVBDivalent(std::istream& istr);
  virtual ~TrialAddAVBDivalent() {}
  //@}
};

inline std::shared_ptr<TrialAddAVBDivalent> MakeTrialAddAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialAddAVBDivalent>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
