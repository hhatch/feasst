"""
Flat-histogram simulation of EMP2 CO2 (https://doi.org/10.1021/j100031a034) in the grand canonical ensemble.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import macrostate_distribution

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/co2_epm2_rigid.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--temperature', type=float, default=298, help='temperature in Kelvin')
PARSER.add_argument('--beta_mu', type=float, default=-6, help='beta time chemical potential')
PARSER.add_argument('--cutoff', type=float, default=12, help='real space cutoff distance')
PARSER.add_argument('--max_particles', type=int, default=300, help='maximum number of particles')
PARSER.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
PARSER.add_argument('--min_particles_second_window', type=int, default=45, help='minimum number of particles in the second window')
PARSER.add_argument('--min_sweeps', type=int, default=2,
                    help='Minimum number of sweeps defined in https://dx.doi.org/10.1063/1.4918557')
PARSER.add_argument('--cubic_side_length', type=float, default=28,
                    help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=0,
                    help='number of iterations for equilibration')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.2, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['prefix'] = 'co2'
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.95*PARAMS['hours_terminate'] - 0.05 # terminate FEASST before SLURM
PARAMS['hours_terminate'] *= PARAMS['procs_per_node'] # real time -> cpu time
PARAMS['hours_checkpoint'] *= PARAMS['procs_per_node']
PARAMS['num_sims'] = PARAMS['num_nodes']
PARAMS['procs_per_sim'] = PARAMS['procs_per_node']
PARAMS['beta'] = 1./(PARAMS['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ
PARAMS['mu'] = PARAMS['beta_mu']/PARAMS['beta']
PARAMS['dccb_cut'] = 4.
PARAMS['dccb_cut'] = PARAMS['cubic_side_length']/int(PARAMS['cubic_side_length']/PARAMS['dccb_cut']) # maximize inside box
PARAMS['ewald_alpha'] = 5.6/PARAMS['cubic_side_length']
PARAMS['mu_init']=10
PARAMS['angle_width'] = 0.01
PARAMS['angle_center'] = np.pi - 0.5*PARAMS['angle_width']

# write TrialGrowFile
with open(PARAMS['prefix']+'_grow_grand_canonical.txt', 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 1 reference_index 0
bond true mobile_site 1 anchor_site 0 num_steps 1 reference_index 0
angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0
""")
with open(PARAMS['prefix']+'_grow_canonical.txt', 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 weight 0.1 angle true mobile_site 2 anchor_site 0 anchor_site2 1 reference_index 0

particle_type 0 weight 0.1 angle true mobile_site 1 anchor_site 0 anchor_site2 2 reference_index 0
""")

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} min0 {min_particles} min1 {min_particles_second_window} num {procs_per_node} overlap 0 alpha 2.15 min_size 3
Checkpoint checkpoint_file {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_side_length {cubic_side_length} particle_type0 {fstprt} cutoff {cutoff} sigma0_1 2.89170901025674 group0 oxygen oxygen_site_type 0
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
RefPotential VisitModel DontVisitModel
#RefPotential Model HardSphere group oxygen cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut} cell_group oxygen
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.25 particle_type 0 tunable_param 0.2 tunable_target_acceptance 0.25
TrialGrowFile grow_file {prefix}_grow_canonical.txt
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialGrowFile grow_file {prefix}_grow_grand_canonical.txt
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.txt include_bonds true
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial name_contains add
RemoveTrial name_contains remove
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} min_flatness 25 collect_flatness 20 min_collect_sweeps 1
TrialGrowFile grow_file {prefix}_grow_grand_canonical.txt
Log trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].txt include_bonds true
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_eq.xyz stop_after_iteration 1
Movie trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index].xyz start_after_iteration 1
Tune trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
#AnalyzeBonds trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_bonds.txt start_after_iteration 1 angle_bin_center {angle_center} angle_bin_width {angle_width} bond_bin_center 1.149 bond_bin_width 0.00001
#PairDistribution trials_per_update 1000 trials_per_write {trials_per_iteration} start_after_iteration 1 \
  dr 0.025 output_file {prefix}n{node}s[sim_index]_gr.txt multistate true multistate_aggregate false
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} output_file {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    lnp = macrostate_distribution.splice_collection_matrix(prefix=params['prefix']+'n0s', suffix='_crit.txt', use_soft=True)
#    lnp.reweight(-0.75, inplace=True)
#    delta_beta_mu = lnp.equilibrium(delta_beta_mu_guess=0.01)
#    #print(delta_beta_mu)
#    #print(lnp.dataframe())
#    #lnp.plot()
#    #plt.savefig(params['prefix']+'.png')
#    #print(lnp.minimums())
#    vapor, liquid = lnp.split()
#    volume = params['cubic_side_length']**3
#    na = physical_constants.AvogadroConstant().value()
#    dens_conv = 1./volume/na*44.01/1e3*1e30 # convert from N/V units of molecules/A^3 to kg/m^3
#    print('vapor density(kg/m^3)', vapor.average_macrostate()*dens_conv)
#    print('liquid density(kg/m^3)', liquid.average_macrostate()*dens_conv)
#    assert np.abs(256 - vapor.average_macrostate()*dens_conv) < 10
#    assert np.abs(712 - liquid.average_macrostate()*dens_conv) < 50

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
