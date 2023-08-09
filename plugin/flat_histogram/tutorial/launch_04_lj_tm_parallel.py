"""
Flat-histogram simulation of single-site Lennard Jones particles in the grand canonical ensemble.
The default temperature is above the critical point.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio

# Parse arguments from command line or change their default values.
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--feasst_install', type=str, default=os.path.expanduser('~')+'/feasst/build/',
    help='FEASST install directory (e.g., the path to build)')
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/lj.fstprt',
    help='FEASST particle definition')
parser.add_argument('--beta', type=float, default=1./1.5, help='inverse temperature')
parser.add_argument('--mu', type=float, default=-2.352321, help='chemical potential')
parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
parser.add_argument('--max_particles', type=int, default=370, help='maximum number of particles')
parser.add_argument('--min_particles', type=int, default=0, help='minimum number of particles')
parser.add_argument('--min_sweeps', type=int, default=1e4, help='minimum number of particles')
parser.add_argument('--cubic_box_length', type=float, default=8,
    help='cubic periodic boundary length')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e6),
    help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibraiton')
parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=5*24, help='hours until termination')
parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
parser.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
parser.add_argument('--run_type', '-r', type=int, default=0,
    help='0: run, 1: submit to queue, 2: post-process')
parser.add_argument('--seed', type=int, default=-1,
    help='Random number generator seed. If -1, assign random seed to each sim.')
parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
parser.add_argument('--scratch', type=str, default=None,
    help='Optionally write scheduled job to scratch/logname/jobid.')
parser.add_argument('--node', type=int, default=0, help='node ID')
parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
args, unknown_args = parser.parse_known_args()
assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
params = vars(args)
params['script'] = __file__
params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
params['hours_terminate'] = 0.99*params['hours_terminate']*params['procs_per_node'] - 0.0333 # terminate before queue
params['hours_checkpoint'] *= params['procs_per_node']
params['num_sims'] = params['num_nodes']
params['procs_per_sim'] = params['procs_per_node']

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 2.25 min_size 5
Checkpoint file_name {prefix}n{node}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.txt
Tune
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialTransfer weight 2 particle_type 0
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].xyz
Tune trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def post_process(params):
    # skip the following checks if temperature is not 1.5
    if np.abs(params['beta'] - 1./1.5) > 1e-5:
        return
    lnpi=pd.read_csv(params['prefix']+'n0_lnpi.txt')
    self.assertAlmostEqual(310.4179421879679, (np.exp(lnpi["ln_prob"]) * lnpi["state"]).sum(), delta=0.5)
    fst = pd.read_csv(params['prefix']+'_lnpi.csv')
    srsw = pd.read_csv('../test/data/stat150.csv')
    plt.plot(fst['state'], fst['ln_prob'], label='FEASST')
    plt.plot(srsw['N'], srsw['lnPI'], linestyle='dashed', label='SRSW')
    plt.xlabel('number of particles', fontsize=16)
    plt.ylabel('ln probability', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(params['prefix']+'_lnpi.png', bbox_inches='tight', transparent='True')

if __name__ == '__main__':
    feasstio.run_simulations(params=params,
                             sim_node_dependent_params=None,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=args)
