"""
This tutorial is similar to tutorial 4, but for this low temperature simulation,
we will split the simulation into two different nodes.
The first node will have less particles but a higher number of sweeps required.
The second node will have dccb but not avb.
"""

import os
import subprocess
import argparse
import random
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
parser.add_argument('--beta', type=float, default=1./0.7, help='inverse temperature')
parser.add_argument('--mu', type=float, default=-4.1603632, help='chemical potential')
parser.add_argument('--mu_init', type=float, default=10, help='initial chemical potential')
parser.add_argument('--num_particles', type=int, default=475, help='number of particles')
parser.add_argument('--num_particles_first_node', type=int, default=375,
  help='number of particles in the first node')
parser.add_argument('--cubic_box_length', type=float, default=8,
    help='cubic periodic boundary length')
parser.add_argument('--dccb_cut', type=float, default=1,
    help='dual-cut configurational bias cutoff')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e6),
    help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibraiton')
parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=5*24, help='hours until termination')
parser.add_argument('--procs_per_node', type=int, default=32, help='number of processors')
parser.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
parser.add_argument('--run_type', '-r', type=int, default=0,
    help='0: run, 1: submit to queue, 2: post-process')
parser.add_argument('--seed', type=int, default=-1,
    help='Random number generator seed. If -1, assign random seed to each sim.')
parser.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
parser.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
parser.add_argument('--node', type=int, default=0, help='node ID')
parser.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
parser.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
# Define sim-dependent parameters in run(sim, ...), with sim integer range of [0, num_sims-1].
args, unknown_args = parser.parse_known_args()
assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
params = vars(args)
params['script'] = __file__
params['sim_id_file'] = params['prefix']+ '_sim_ids.txt'
params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
params['num_sims'] = params['num_nodes']
params['dccb_cut'] = params['cubic_box_length']/int(params['cubic_box_length']/params['dccb_cut'])

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_checkpoint} ln_prob_file {prefix}n{node}_lnpi.txt min_window_size -1
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint file_name {prefix}n{node}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
NeighborCriteria maximum_distance 1.375 minimum_distance 0.9
{lj_potential}
{ref_potential}
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
{avb_trials}
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization and nvt equilibration
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_eq.txt
Tune
TrialAdd particle_type 0
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
{gce_trial}
Log trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index].xyz
Tune trials_per_write {trials_per_iteration} file_name {prefix}en{node}s[sim_index]_tune.txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_en.txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}n{node}s[sim_index]_crit.txt
""".format(**params))

def run(sim, params):
    """ Run a single simulation. If all simulations are complete, run PostProcess. """
    if args.queue_task == 0:
        params['sim'] = sim + params['node']*params['procs_per_node']
        if params['seed'] == -1:
            params['seed'] = random.randrange(int(1e9))
        if params['node'] == 0:
            params['min_particles']=0
            params['max_particles']=params['num_particles_first_node']
            params['gce_trial']='TrialTransfer weight 2 particle_type 0\nTrialTransferAVB weight 0.2 particle_type 0'
            params['lj_potential']='Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJones'
            params['ref_potential']=''
            params['avb_trials']='TrialAVB2 weight 0.1 particle_type 0\nTrialAVB4 weight 0.1 particle_type 0'
            params['min_sweeps']=2000
            params['window_alpha']=2
            params['min_window_size']=5
        elif params['node'] == 1:
            params['min_particles']=params['num_particles_first_node']
            params['max_particles']=params['num_particles']
            params['gce_trial'] = 'TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 10'
            params['lj_potential']='Potential Model LennardJones'
            params['ref_potential']="""RefPotential Model LennardJones cutoff {dccb_cut} VisitModel VisitModelCell min_length {dccb_cut}""".format(**params)
            params['avb_trials']=''
            params['min_sweeps']=200
            params['window_alpha']=1
            params['min_window_size']=3
        file_name = params['prefix']+str(sim)+'_launch_run'
        write_feasst_script(params, file_name=file_name+'.txt')
        syscode = subprocess.call(
            args.feasst_install+'bin/fst < '+file_name+'.txt  > '+file_name+'.log',
            shell=True, executable='/bin/bash')
    else: # if queue_task < 1, restart from checkpoint
        syscode = subprocess.call(
            args.feasst_install+'bin/rst '+params['prefix']+str(sim)+'_checkpoint.fst',
            shell=True, executable='/bin/bash')
    if syscode == 0: # if simulation finishes with no errors, write to sim id file
        with open(params['sim_id_file'], 'a', encoding='utf-8') as file1:
            file1.write(str(sim)+'\n')
        # if all sims are complete, post process or test once (by removing sim id file)
        if feasstio.all_sims_complete(params['sim_id_file'], params['num_sims']):
            os.remove(params['sim_id_file'])
            post_process(params)
    return syscode

def post_process(params):
    assert False

if __name__ == '__main__':
    feasstio.run_simulations(params=params,
                             run_function=run,
                             post_process_function=post_process,
                             queue_function=feasstio.slurm_single_node,
                             run_type=args.run_type,
                             queue_id=args.queue_id,
                             queue_task=args.queue_task)
