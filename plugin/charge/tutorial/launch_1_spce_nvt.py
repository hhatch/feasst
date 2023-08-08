"""
Example SPC/E canonical ensemble Monte Carlo simulation using FEASST.
Approximately compare with https://doi.org/10.1063/1.476834.
There are systematic differences in the energy due to different Ewald cutoffs, etc.
"""

import os
import subprocess
import argparse
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import feasstio
from pyfeasst import physical_constants

# Parse arguments from command line or change their default values.
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--feasst_install', type=str, default=os.path.expanduser('~')+'/feasst/build/',
    help='FEASST install directory (e.g., the path to build)')
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/spce.fstprt',
    help='FEASST particle definition')
parser.add_argument('--temperature', type=float, default=298, help='temperature in Kelvin')
parser.add_argument('--num_particles', type=int, default=512, help='number of particles')
parser.add_argument('--cubic_box_length', type=float, default=24.8586887,
    help='cubic periodic boundary length')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e5),
    help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibraiton')
parser.add_argument('--production_iterations', type=int, default=int(1e1),
    help='number of iterations for production')
parser.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
parser.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
parser.add_argument('--prefix', type=str, default='spce', help='prefix for all output file names')
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
params['num_sims'] = params['num_nodes']*params['procs_per_node']
params['cutoff'] = 0.5*params['cubic_box_length']
params['alpha'] = 5.6/params['cubic_box_length']
params['beta'] = 1./(params['temperature']*physical_constants.MolarGasConstant().value()/1e3) # mol/kJ


def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} physical_constants CODATA2010 cutoff {cutoff}
Potential VisitModel Ewald alpha {alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened VisitModel VisitModelCutoffOuter table_size 1e6
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential 1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.txt
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log

# canonical ensemble production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}.xyz
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_en.txt
CPUTime trials_per_write {trials_per_iteration} file_name {prefix}{sim}_cpu.txt
Run until_criteria_complete true
""".format(**params))

def run(sim, params):
    """ Run a single simulation. If all simulations are complete, run PostProcess. """
    if args.queue_task == 0:
        params['sim'] = sim + params['node']*params['procs_per_node']
        if params['seed'] == -1:
            params['seed'] = random.randrange(int(1e9))
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
    """ Approximately compare energy with https://doi.org/10.1063/1.476834 """
    if params['num_sims'] == 1: # compare energy
      log = pd.read_csv(params['prefix']+'0.txt')
      assert int(log['num_particles_of_type0'][0]) == params['num_particles']
      energy = pd.read_csv(params['prefix']+'0_en.txt')
      diff = energy['average'][0] - (-46.82)*params['num_particles']
      assert np.abs(diff) < 20*np.sqrt(energy['block_stdev'][0]**2 + (0.02*params['num_particles'])**2)

if __name__ == '__main__':
    feasstio.run_simulations(params=params,
                             run_function=run,
                             post_process_function=post_process,
                             queue_function=feasstio.slurm_single_node,
                             run_type=args.run_type,
                             queue_id=args.queue_id,
                             queue_task=args.queue_task)
