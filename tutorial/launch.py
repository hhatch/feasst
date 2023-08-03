"""
Example single-site Lennard-Jones canonical ensemble Monte Carlo simulation using FEASST.
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

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cubic_box_length', type=float, default=8, help='cubic periodic boundary length')
parser.add_argument('--fstprt', type=str, default='/feasst/forcefield/lj.fstprt', help='FEASST particle definition')
parser.add_argument('--beta_lower', type=float, default=0.8, help='lowest inverse temperature')
parser.add_argument('--beta_upper', type=float, default=1.2, help='highest inverse temperature')
parser.add_argument('--num_particles', type=int, default=50, help='number of particles')
parser.add_argument('--trials_per_iteration', type=int, default=int(1e5), help='like cycles, but not necessary num_particles')
parser.add_argument('--equilibration_iterations', type=int, default=int(1e1))
parser.add_argument('--production_iterations', type=int, default=int(1e3))
parser.add_argument('--seed', type=int, default=random.randrange(int(1e9)), help='random number generator seed')
parser.add_argument('--hours_checkpoint', type=float, default=1., help='number of hours per checkpoint')
parser.add_argument('--hours_terminate', type=float, default=1., help='number of hours until termination')
parser.add_argument('--sim', type=int, default=0, help='simulation id')
parser.add_argument('--num_procs', type=int, default=2, help='number of processors')
parser.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
parser.add_argument('--run_type', '-r', type=int, default=0, help='0: submit to HPC queue, 1: run, 2: post-process')
parser.add_argument('--slurm_id', type=int, default=-1, help='Automatically input by slurm scheduler. If != -1, read args from file.')
parser.add_argument('--slurm_task', type=int, default=0, help='Automatically input by slurm scheduler. If > 0, restart from checkpoint.')
args = parser.parse_args()
params = vars(args)
params['script'] = __file__ # store for reference
params['minutes'] = round(params['hours_terminate']*60) # incase decimal number of hours
params['hours_checkpoint'] *= params['num_procs'] # real time to CPU time
params['hours_terminate'] = 0.95*params['hours_terminate']*params['num_procs'] # terminate 95% early so SLURM doesn't terminate first
print('params', params)

# write fst script to run a single simulation
def mc(params, file_name):
    with open(file_name, 'w') as myfile: myfile.write("""
# high temperature gcmc to generate initial configuration
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2

# grand canonical ensemble initalization
TrialAdd particle_type 0
Checkpoint file_name {prefix}_checkpoint{sim}.fst num_hours 0.1 num_hours_terminate {hours_terminate}
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
Run until_criteria_complete true
""".format(**params))

# run a single simulation
def run(sim):
    if args.slurm_task == 0:
        params['sim'] = sim
        params['beta'] = np.linspace(params['beta_lower'], params['beta_upper'], num=params['num_procs'])[sim]
        params['seed'] = random.randrange(int(1e9))
        file_name = params['prefix']+str(sim)+'_launch_run'
        mc(params, file_name=file_name+'.txt')
        syscode = subprocess.call('../build/bin/fst < ' + file_name + '.txt  > ' + file_name + str(sim)+'.log', shell=True, executable='/bin/bash')
    else: # if slurm_task < 1, restart from checkpoint
        syscode = subprocess.call('../build/bin/rst checkpoint' + str(sim) + '.fst', shell=True, executable='/bin/bash')
    if sim == 0 and syscode == 0: # if simulation finishes with no errors, post-process
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

# after the simulation is complete, perform some tests or analysis
class PostProcess(unittest.TestCase):
    def test(self):
        import pandas as pd
        log = pd.read_csv('lj0.txt')
        print(log)
        self.assertTrue(log['num_particles_of_type0'][0] == params['num_particles'])

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open(params['prefix'] + '_slurm.txt', 'w') as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {num_procs} -N 1 -t {minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python {script} --run_type 1 --slurm_id $SLURM_JOB_ID --slurm_task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

if __name__ == '__main__':
    if args.run_type == 0: # queue on SLURM
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 " + params['prefix'] + "_slurm.txt | awk '{print $4}' >> " + params['prefix']+ "_launch_ids.txt", shell=True, executable='/bin/bash')
        with open(params['prefix']+ '_launch_ids.txt') as file1:
            slurm_id = file1.readlines()[-1]
        with open('lj_params'+slurm_id, 'w') as file1:
            file1.write(json.dumps(params))
    elif args.run_type == 1: # run directly
        if args.slurm_id != -1: # read parameter from file if submit via slurm
            with open('lj_params'+slurm_id, 'r') as file1:
                params = json.load(file1)
        with Pool(params['num_procs']) as pool:
            codes = pool.starmap(run, zip(range(0, params['num_procs'])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif args.run_type == 2: # post process
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
