"""
Example single-site Lennard-Jones canonical ensemble Monte Carlo simulation using FEASST.
Run multiple densities using multiple processors/nodes/restarts, and plot results.
Compare with T*=0.9 in https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.
Usage: python /path/to/feasst/tutorial/launch.py --help
"""

import sys
import subprocess
import argparse
import unittest
import random
import json
from multiprocessing import Pool
import numpy as np
import pandas as pd
from pyfeasst import feasstio

# parse arguments from command line or change their default values.
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/lj.fstprt', help='FEASST particle definition')
parser.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
parser.add_argument('--num_particles', type=int, default=500, help='number of particles')
parser.add_argument('--rho_lower', type=float, default=1e-3, help='lowest number density')
parser.add_argument('--rho_upper', type=float, default=9e-3, help='highest number density')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e5), help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1), help='number of iterations for equilibraiton')
parser.add_argument('--production_iterations', type=int, default=int(1e3), help='number of iterations for production')
parser.add_argument('--seed', type=int, default=random.randrange(int(1e9)), help='random number generator seed')
parser.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=0.1, help='number of hours until termination')
parser.add_argument('--sim', type=int, default=0, help='simulation ID')
parser.add_argument('--num_procs', type=int, default=5, help='number of processors')
parser.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
parser.add_argument('--run_type', '-r', type=int, default=0, help='0: run, 1: submit to queue, 2: post-process')
parser.add_argument('--slurm_id', type=int, default=-1, help='Automatically input by slurm scheduler. If != -1, read args from file')
parser.add_argument('--slurm_task', type=int, default=0, help='Automatically input by slurm scheduler. If > 0, restart from checkpoint')
parser.add_argument('--max_restarts', type=int, default=10, help='Number of SLURM restarts')
parser.add_argument('--num_nodes', type=int, default=1, help='Number of SLURM nodes')
parser.add_argument('--node', type=int, default=0, help='node ID')
args = parser.parse_args()
params = vars(args)
params['script'] = __file__
params['minutes'] = int(params['hours_terminate']*60) # minutes used for SLURM queue time
params['hours_terminate'] = 0.99*params['hours_terminate'] - 0.0333 # terminate before SLURM
params['num_sims'] = params['num_nodes']*params['num_procs']
params['rhos'] = np.linspace(params['rho_lower'], params['rho_upper'], num=params['num_sims'])
params['cubic_box_lengths'] = np.power(params['num_particles']/params['rhos'], 1./3.).tolist()
params['rhos'] = params['rhos'].tolist()

# write fst script for a single simulation with params keys {} enclosed
def mc(params, file_name):
    with open(file_name, 'w') as myfile: myfile.write("""
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

# run a single simulation
def run(sim):
    if args.slurm_task == 0:
        params['sim'] = sim + params['node']*params['num_procs']
        params['cubic_box_length'] = params['cubic_box_lengths'][sim]
        params['seed'] = random.randrange(int(1e9))
        file_name = params['prefix']+str(sim)+'_launch_run'
        mc(params, file_name=file_name+'.txt')
        syscode = subprocess.call('../build/bin/fst < ' + file_name + '.txt  > ' + file_name + '.log', shell=True, executable='/bin/bash')
    else: # if slurm_task < 1, restart from checkpoint
        syscode = subprocess.call('../build/bin/rst ' + params['prefix'] + str(sim) + '_checkpoint.fst', shell=True, executable='/bin/bash')
    if syscode == 0: # if simulation finishes with no errors, write to sim id file
        with open(params['sim_id_file'], 'a') as file1:
            file1.write(str(sim)+'\n')
        # post process / test if all sims are complete (ensure once by clearing sim id file)
        if feasstio.all_sims_complete(params['sim_id_file'], params['num_sims']):
            open(params['sim_id_file'], 'w').close() # clear file
            unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

# after the simulation is complete, compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
class PostProcess(unittest.TestCase):
    def test(self):
        ens = np.zeros(shape=(params['num_sims'], 2))
        for sim in range(params['num_sims']):
            log = pd.read_csv('lj'+str(sim)+'.txt')
            self.assertTrue(int(log['num_particles_of_type0'][0]) == params['num_particles'])
            en = pd.read_csv('lj'+str(sim)+'_en.txt')
            ens[sim] = np.array([en['average'][0], en['block_stdev'][0]])/params['num_particles']
        import matplotlib.pyplot as plt
        # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        rhos_srsw = [0.001, 0.003, 0.005, 0.007, 0.009]
        ens_srsw = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02]
        en_stds_srsw = [1.89E-05, 3.21E-05, 3.80E-05, 7.66E-05, 2.44E-05]
        plt.errorbar(rhos_srsw, ens_srsw, en_stds_srsw, fmt='+', label='SRSW')
        plt.errorbar(params['rhos'], ens[:, 0], ens[:, 1], fmt='x', label='FEASST')
        plt.xlabel(r'$\rho$', fontsize=16)
        plt.ylabel(r'$U/N$', fontsize=16)
        plt.legend(fontsize=16)
        plt.savefig(params['prefix']+'_en.png', bbox_inches='tight', transparent='True')
        if len(rhos_srsw) == params['num_sims']: # compare with srsw exactly
            for sim in range(params['num_sims']):
                self.assertAlmostEqual(ens[sim][0], ens_srsw[sim],
                    delta=1.96*np.sqrt(ens[sim][1]**2 + en_stds_srsw[sim]**2))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open(params['prefix'] + '_slurm.txt', 'w') as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {num_procs} -N 1 -t {minutes}:00 -o {prefix}_slurm_%j.txt -e {prefix}_slurm_%j.txt
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python {script} --run_type 0 --node {node} --slurm_id $SLURM_ARRAY_JOB_ID --slurm_task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

if __name__ == '__main__':
    params['sim_id_file'] = params['prefix']+ "_sim_ids.txt"
    open(params['sim_id_file'], 'w').close() # clear file, then append sim id when complete
    if args.run_type == 0: # run directly
        if args.slurm_id != -1: # if run from SLURM
            if args.slurm_task == 0: # read param file if not checkpoint
                with open('lj_params'+str(args.slurm_id)+'.json', 'r') as file1:
                    params = json.load(file1)
        else:
            with open('lj_params.json', 'w') as file1:
                file1.write(json.dumps(params))
        with Pool(params['num_sims']) as pool:
            codes = pool.starmap(run, zip(range(0, params['num_sims'])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif args.run_type == 1: # queue on SLURM
        params['slurm_id_file'] = params['prefix']+ "_slurm_ids.txt"
        open(params['slurm_id_file'], 'w').close() # empty file contents
        for node in range(params['num_nodes']):
            params['node'] = node
            slurm_queue()
            subprocess.call("sbatch --array=0-" + str(params['max_restarts']) + "%1 " + params['prefix'] + "_slurm.txt | awk '{print $4}' >> " + params['slurm_id_file'], shell=True, executable='/bin/bash')
            with open(params['slurm_id_file'], 'r') as file1:
                slurm_id = file1.read().splitlines()[-1]
            with open('lj_params'+slurm_id+'.json', 'w') as file1:
                file1.write(json.dumps(params))
    elif args.run_type == 2: # post process
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
