import pandas as pd
import feasst as fst
import pyfeasst

# read checkpoints
num_procs = 4
clones = fst.MakeClones('checkpoint', num_procs)

# create a temporary aggregate criteria for reweighting
#criteria = fst.FlatHistogram(
#    fst.MacrostateNumParticles(fst.Histogram(fst.args({"width": "1",
#        "min": str(int(clones.flat_histogram(0).macrostate().histogram().center_of_bin(0))),
#        "max": str(clones.flat_histogram(num_procs-1).macrostate().histogram().center_of_last_bin())}))),
#    fst.MakeTransitionMatrix(fst.args({"min_sweeps": "0"})),
#    fst.args({"beta": str(clones.clone(0).criteria().beta()),
#              "chemical_potential0": str(clones.clone(0).criteria().chemical_potential(0))}))
#criteria.set_ln_prob(clones.ln_prob())
gce = fst.GrandCanonicalEnsemble(clones)

# reweight criteria
#dbetamurw = -1
#gce = fst.GrandCanonicalEnsemble(criteria)
gce.reweight(-1)
#lnpi_rw = criteria.reweight(dbetamurw)
#criteria.set_ln_prob(lnpi_rw)
#criteria.set_chemical_potential(criteria.chemical_potential()+dbetamurw/criteria.beta())

import matplotlib.pyplot as plt
plt.plot(gce.ln_prob().values())
plt.title(r'$\beta\mu=$'+ str(gce.beta_mu()))
#plt.show()

# read the energy moments
#num_analyzers = clones.clone(0).num_analyzers()
#num_index = num_analyzers - 2
#assert(clones.clone(0).analyze(num_index).class_name() == "AnalyzeFactory")
#assert(clones.clone(0).analyze(num_index).analyze(0).class_name() == "Energy")
#print(type(u_moments))
#u_moments.resize(2, fst.DoubleVector())
#u_moments[0] = fst.DoubleVector()
#print(type(u_moments[0]))
#u_moments[0].resize(1)
#u_moments[1] = fst.DoubleVector()
#u_moments[1].resize(1)
u_nvt = fst.DoubleVector()
clones.stitch(u_nvt, "Energy", fst.AccumulatorMoment(0));
#clones.stitch(u_moments[0], "Energy", fst.AccumulatorMoment(0));
u2_nvt = fst.DoubleVector()
clones.stitch(u2_nvt, "Energy", fst.AccumulatorMoment(1));
#clones.stitch(u_moments[1], "Energy", fst.AccumulatorMoment(1));
u_moments = fst.Double2DVector(2)
u_moments[0] = u_nvt
u_moments[1] = u2_nvt
print(type(u_moments[0]))
#nu = list()
#for num, energy in enumerate(u_nvt): nu.append(num*energy)
#
#gc_u_nvt = gce.average(u_nvt)
#gc_u2_nvt = gce.average(u2_nvt)
#gc_nu = gce.average(nu)
#gc_n = gce.average_macrostate()

# extrapolate
#ln_prob_new = clones.ln_prob()
#u_nvt_new = u_nvt
#mu = clones.clone(0).criteria().chemical_potential(0)
#dbeta = 1./1.2 - 1./1.5
#for state, lnp in enumerate(ln_prob_new.values()):
#    dlnpdb = -u_nvt[state] + gc_u_nvt + mu*state
#    d2lnpdb2 = u2_nvt[state] - u_nvt[state]**2 - gc_u2_nvt + gc_u_nvt**2 - mu*(gc_nu-gc_n*gc_u_nvt)
#    ln_prob_new.add(state,  dlnpdb*dbeta + d2lnpdb2*dbeta**2/2.)
#    dunvtdb = -u2_nvt[state] + u_nvt[state]**2
#    u_nvt_new[state] += dunvtdb*dbeta
#    #<u^2> - <u>^2 - gc<u^2> + gc<u>^2 - mu*(gc<nu>-gc<n>gc<u>)

gce.extrapolate_beta(u_moments,
    fst.args({"beta_new": "1.2",
              "beta_original": str(clones.clone(0).criteria().beta())}));
#ln_prob_new.normalize()
plt.plot(gce.ln_prob().values())
#plt.show()

#criteria.set_ln_prob(ln_prob_new)
#criteria.set_beta(criteria.beta() + dbeta)
#gce = fst.GrandCanonicalEnsemble(criteria)
print('betamu', gce.beta_mu())

# find saturation at extrapolated conditions
#gce.set_ln_prob(ln_prob_new)
gce=pyfeasst.find_equilibrium(gce)

volume = clones.clone(0).configuration().domain().volume()
num_vapor = gce.average_macrostate(0)
num_liquid = gce.average_macrostate(1)
print('vapor density', num_vapor/volume)
print('liquid density', num_liquid/volume)
print('pressure', gce.betaPV()/volume/clones.clone(0).criteria().beta())
print('vapor en', gce.average(u_moments[0], 0)/num_vapor)
print('liquid en', gce.average(u_moments[0], 1)/num_liquid)
#print('betamu', sat.beta_mu())

plt.plot(gce.ln_prob().values())
plt.show()

