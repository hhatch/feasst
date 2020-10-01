import copy
import matplotlib.pyplot as plt
import feasst as fst
import pyfeasst

for temp in [1.2, 1, 0.8]:
    # read checkpoints
    num_procs = 4
    clones = fst.MakeClones('checkpoint', num_procs)
    gce = fst.GrandCanonicalEnsemble(clones)
    if temp == 1.2: plt.plot(gce.ln_prob().values(),
                             label='T*='+str(1./clones.clone(0).criteria().beta()))

    # collect the energy moments
    u_moments = fst.Double2DVector(2)
    for moment in range(2):
        u = fst.DoubleVector()
        clones.stitch(u, "Energy", fst.AccumulatorMoment(moment))
        u_moments[moment] = u

    gce.extrapolate_beta(u_moments,
        fst.args({"beta_new": str(1/temp),
                  "beta_original": str(clones.clone(0).criteria().beta())}))
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

    plt.plot(gce.ln_prob().values(), label='T*='+str(temp))

plt.legend()
plt.show()
