import pandas as pd
import feasst as fst
import pyfeasst

# read checkpoints
num_procs = 12
clones = fst.MakeClones('checkpoint', num_procs)
volume = clones.clone(0).configuration().domain().volume()
beta = clones.clone(0).criteria().beta()

gce = pyfeasst.find_equilibrium(fst.GrandCanonicalEnsemble(clones))
#plt.plot(sat.bias().ln_prob().values())
#plt.title(r'$\beta\mu=$'+ str(sat.beta_mu(0)))
#plt.show()

# compute saturation pressure
R=clones.clone(0).configuration().physical_constants().ideal_gas_constant()
na=clones.clone(0).configuration().physical_constants().avogadro_constant()
press_conv=R/1e3*1e30/na
print('saturation pressure (kPa)', gce.betaPV()/volume/beta*press_conv)

# compute saturation compositions
num_index = fst.SeekAnalyze().index("NumParticles", clones.clone(0))[0]
num0 = list()
for iclone in range(clones.num()):
    num_states = clones.flat_histogram(iclone).num_states()
    for state in range(num_states):
        if iclone == clones.num() - 1 or state < num_states - 3:  # ignore last 3 states for all but last processor
            num0.append(clones.clone(iclone).analyze(num_index).analyze(state).accumulator().average())
num_vapor = gce.average(num0, 0)
num_liquid = gce.average(num0, 1)
print(fst.SeekAnalyze().multistate_average("NumParticles", clones.clone(0)))
print(len(fst.SeekAnalyze().multistate_average("NumParticles", clones.clone(0))))

print('vapor y_C2H4', 1 - num_vapor/gce.average_macrostate(0))
print('liquid x_C2H4', 1 - num_liquid/gce.average_macrostate(1))

# obtain extensive moments
extmom_index = fst.SeekAnalyze().index("ExtensiveMoments", clones.clone(0))[0]

# Return ExtensiveMoments, dervied class of Analyze, by serialization, which is relatively slow.
# see steppers/include/ExtensiveMoments.h for moments API
def extensive_moment(window, state):
    return fst.ExtensiveMoments(clones.clone(window).analyze(extmom_index).analyze(state))

extmom = extensive_moment(2, 30)
for p in range(3):
    for m in range(2):
        for k in range(2):
            for j in range(1):
                for i in range(1):
                    print(p, m, k, j, i,
                          extmom.moments(p, m, k, j, i).num_values(),
                          extmom.moments(p, m, k, j, i).sum_dble(),
                          extmom.moments(p, m, k, j, i).sum_of_squared_dble())

