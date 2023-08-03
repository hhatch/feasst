"""
Usage:

# fill an HPC SLURM node with 32 NVT MC LJ simulations at various temperatures, restarting upon the hour.
python /path/to/feasst/tutorial/launch.py

# run manually in your terminal
python /path/to/feasst/tutorial/launch.py -r 1

# run post processing and tests
python /path/to/feasst/tutorial/launch.py -r 2
"""

#HWH: recast per 'cycle' instead of per trial

import sys
import subprocess
import argparse
import unittest
import random
from multiprocessing import Pool

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/lj.fstprt", "beta": 1.2,
    "num_particles": 50,
    "trials_per_iteration": 1e5,
    "equilibration_iterations": 1e6,
    "production_iterations": 1e8,
    "seed": random.randrange(int(1e9)), "num_hours": 1, "script": __file__}
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]

# write fst script to run a single simulation
def mc_lj(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# high temperature gcmc to generate initial configuration
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
TrialAdd particle_type 0
Checkpoint file_name checkpoint{sim}.fst num_hours 0.1 num_hours_terminate {num_hours_terminate}
Run until_num_particles {num_particles}

# nvt equilibration
RemoveTrial name TrialAdd
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
Tune
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-8
Run until_criteria_complete true

# nvt production
RemoveModify name Tune
Log trials_per_write {trials_per_iteration} file_name lj{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name lj{sim}.xyz
Energy trials_per_write {trials_per_iteration} file_name lj_en{sim}.txt
Run num_trials {production_iterations}
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python {script} --run_type 1 --task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

# add additional arguments for multi-core simulations
params.update({"sim": 0, "num_nodes": 1, "procs_per_node": 2})
params["num_sims"] = params["num_nodes"]*params["procs_per_node"]

# set a simulation parameter to vary for each processor
betas = np.linspace(0.8, 1.2, num=params["num_sims"])

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: submit batch to scheduler, 1: run batch on host")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# run a single simulation as part of the batch to fill a node
def run(sim):
    if args.task == 0:
        params["sim"] = sim
        params["beta"] = betas[sim]
        params["seed"] = random.randrange(int(1e9))
        file_name = "launch_run"+str(sim)+".txt"
        mc_lj(params, file_name=file_name)
        syscode = subprocess.call("../build/bin/fst < " + file_name + " > launch_run"+str(sim)+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../build/bin/rst checkpoint" + str(sim) + ".fst", shell=True, executable='/bin/bash')
    if syscode == 0: # if simulation finishes with no errors, post-process
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

# after the simulation is complete, perform some tests or analysis
class PostProcess(unittest.TestCase):
    def test(self):
        import pandas as pd
        log = pd.read_csv('lj0.txt')
        print(df)
        self.assertTrue(log['attempt'][-1] > 1)

if __name__ == "__main__":
    if args.run_type == 0:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(params["num_sims"]) as pool:
            codes = pool.starmap(run, zip(range(0, params["num_sims"])))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
