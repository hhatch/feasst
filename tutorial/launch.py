"""
Example single-site Lennard-Jones canonical ensemble Monte Carlo simulation using FEASST.
Run multiple densities using multiple processors/nodes/restarts, and plot results.
Compare with T*=0.9 in https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.
Usage: python /path/to/feasst/tutorial/launch.py --help
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
parser.add_argument('--feasst_install', type=str, default=os.path.expanduser("~")+'/feasst/build/',
    help='FEASST install directory (e.g., the path to build)')
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/lj.fstprt',
    help='FEASST particle definition')
parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
parser.add_argument('--rho_lower', type=float, default=1e-3, help='lowest number density')
parser.add_argument('--rho_upper', type=float, default=9e-3, help='highest number density')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e5),
    help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1),
    help='number of iterations for equilibraiton')
parser.add_argument('--production_iterations', type=int, default=int(1e3),
    help='number of iterations for production')
parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
parser.add_argument('--num_procs', type=int, default=5, help='number of processors')
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
params['sim_id_file'] = params['prefix']+ "_sim_ids.txt"
params['minutes'] = int(params['hours_terminate']*60) # minutes allocated on queue
params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before queue
params['num_sims'] = params['num_nodes']*params['num_procs']
params['rhos'] = np.linspace(params['rho_lower'], params['rho_upper'], num=params['num_sims'])
params['cubic_box_lengths'] = np.power(params['num_particles']/params['rhos'], 1./3.).tolist()
params['rhos'] = params['rhos'].tolist()

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
# high temperature gcmc to generate initial configuration
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model LennardJones VisitModel VisitModelCell min_length max_cutoff
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours 0.1 num_hours_terminate {hours_terminate}

# grand canonical ensemble initalization
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd

# canonical ensemble equilibration
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Run until_criteria_complete true
RemoveModify name Tune

# nvt production
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
        params['sim'] = sim + params['node']*params['num_procs']
        params['cubic_box_length'] = params['cubic_box_lengths'][sim]
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
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    ens = np.zeros(shape=(params['num_sims'], 2))
    for sim in range(params['num_sims']):
        log = pd.read_csv('lj'+str(sim)+'.txt')
        assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv('lj'+str(sim)+'_en.txt')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
    # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
    rhos_srsw = [0.001, 0.003, 0.005, 0.007, 0.009]
    ens_srsw = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02]
    en_stds_srsw = [1.89E-05, 3.21E-05, 3.80E-05, 7.66E-05, 2.44E-05]
    plt.errorbar(rhos_srsw, ens_srsw, en_stds_srsw, fmt='+', label='SRSW')
    plt.errorbar(params['rhos'], ens[:, 0], ens[:, 1], fmt='x', label='FEASST')
    plt.xlabel(r'$\rho$', fontsize=16)
    plt.ylabel(r'$U/N$', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(params['prefix']+'_energy.png', bbox_inches='tight', transparent='True')
    if len(rhos_srsw) == params['num_sims']: # compare with srsw exactly
        for sim in range(params['num_sims']):
            diff = ens[sim][0] - ens_srsw[sim]
            assert np.abs(diff) < 1.96*np.sqrt(ens[sim][1]**2 + en_stds_srsw[sim]**2)

if __name__ == '__main__':
    feasstio.run_simulations(params=params,
                             run_function=run,
                             post_process_function=post_process,
                             queue_function=feasstio.slurm_single_node,
                             run_type=args.run_type,
                             queue_id=args.queue_id,
                             queue_task=args.queue_task)
