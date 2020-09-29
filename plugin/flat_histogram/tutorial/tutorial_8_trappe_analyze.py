import pandas as pd
import feasst as fst
import pyfeasst

# read checkpoints
clones = fst.MakeClones('checkpoint', 12)

#import matplotlib.pyplot as plt
#plt.plot(clones.ln_prob().values())
#plt.title(r'$\beta\mu=$'+ str(clones.clone(0).criteria().beta_mu(0)))
#plt.show()

# create an aggregate criteria
criteria = fst.FlatHistogram(
    fst.MacrostateNumParticles(fst.Histogram(fst.args({"width": "1", "min": "0", "max": "400"}))),
    fst.MakeTransitionMatrix(fst.args({"min_sweeps": "20"})),
    fst.args({"beta": str(clones.clone(0).criteria().beta()),
              "chemical_potential0": str(clones.clone(0).criteria().chemical_potential(0)),
              "chemical_potential1": str(clones.clone(0).criteria().chemical_potential(1))}))
criteria.set_ln_prob(clones.ln_prob())

# find saturation
def objective(criteria, beta_mu_rw):
    delta_conjugate = beta_mu_rw - criteria.beta_mu()
    ln_prob_rw = criteria.reweight(delta_conjugate)
    return ln_prob_rw.saturation_objective(delta_conjugate)

from scipy.optimize import minimize
def find_sat(criteria, beta_mu_guess=-1):
    res = minimize(lambda beta_mu_rw: objective(criteria, beta_mu_rw[0]), beta_mu_guess, tol=1e-8)
    mu_saturation = res["x"][-1]/criteria.beta()
    delta_conjugate = criteria.beta()*mu_saturation - criteria.beta_mu()
    criteria.set_chemical_potential(mu_saturation)
    criteria.set_ln_prob(criteria.reweight(delta_conjugate))
    return criteria

sat=find_sat(criteria)
#plt.plot(sat.bias().ln_prob().values())
#plt.title(r'$\beta\mu=$'+ str(sat.beta_mu(0)))
#plt.show()

# print saturation pressure
R=clones.clone(0).configuration().physical_constants().ideal_gas_constant()
na=clones.clone(0).configuration().physical_constants().avogadro_constant()
press_conv=R/1e3*1e30/na
print('saturation pressure (kPa)',
      sat.pressure(clones.clone(0).configuration().domain().volume(), 0)*press_conv)

# print saturation compositions
num_analyzers = clones.clone(0).num_analyzers()
num_index = num_analyzers - 4
assert(clones.clone(0).analyze(num_index).class_name() == "AnalyzeFactory")
print(clones.clone(0).analyze(num_index).analyze(3).accumulator().average())

num0 = list()
for iclone in range(clones.num()):
    num_states = clones.flat_histogram(iclone).num_states()
    for state in range(num_states):
        if iclone == clones.num() - 1 or state < num_states - 3:  # ignore last 3 states for processors 0-10
            num0.append(clones.clone(iclone).analyze(num_index).analyze(state).accumulator().average())
num_vapor = sat.average(num0, 0)
num_liquid = sat.average(num0, 1)
#print(num_vapor, num_liquid)

print('vapor y_C2H4', 1 - num_vapor/sat.average_macrostate(0))
print('liquid x_C2H4', 1 - num_liquid/sat.average_macrostate(1))

# obtain extensive moments
extmom_index = num_analyzers - 3
assert(clones.clone(0).analyze(extmom_index).class_name() == "ExtensiveMoments")
# HWH fix this: one for each macrostate
#extmom = fst.ExtensiveMoments(clones.clone(0).analyze(extmom_index))
#print(extmom.moments(2, 0, 0, 0, 0).sum_dble())
#print(extmom.moments(2, 0, 0, 0, 0).sum_of_squared_dble())

