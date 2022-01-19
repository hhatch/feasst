import subprocess
import numpy as np
import argparse
from multiprocessing import Pool
import random

# define parameters of a pure component NVT MC Lennard-Jones simulation
params = {
    "length": 8, "fstprt": "/feasst/forcefield/lj.fstprt", "beta": 1.2,
    "num_particles": 350, "equilibration": 1e6, "production": 1e8, # slow
    #"num_particles": 50, "equilibration": 1e6, "production": 1e6, # fast
    "steps_per": 1e5, "seed": "time",
    "sim": 0, "num_hours": 0.05, "num_nodes": 1, "procs_per_node": 3}
params["num_sims"] = params["num_nodes"]*params["procs_per_node"]
params["num_minutes"] = round(params["num_hours"]*60)
params["num_hours_terminate"] = 0.95*params["num_hours"]

# write fst script to run a single simulation
def mc_lj(params=params, file_name="tutorial.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# high temperature gcmc to generate initial configuration
RandomMT19937 seed {seed}
Configuration cubic_box_length {length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta 0.1 chemical_potential 10
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
TrialAdd particle_type 0
Checkpoint file_name checkpoint{sim}.fst num_hours 0.1 num_hours_terminate {num_hours_terminate}
Run until_num_particles {num_particles}

# nvt equilibration
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Tune steps_per {steps_per}
CheckEnergy steps_per {steps_per} tolerance 1e-8
Run num_attempts {equilibration}

# nvt production
RemoveModify name Tune
Log steps_per {steps_per} file_name lj{sim}.txt
Movie steps_per {steps_per} file_name lj{sim}.xyz
Energy steps_per_write {steps_per} file_name en{sim}.txt
Run num_attempts {production}
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
python tutorial.py --run_type 1 --task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

# set a simulation parameter to vary for each processor
betas = np.linspace(0.8, 1.2, num=params["num_sims"])

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: run single simulation on host, 1: run batch on host, 2: submit match to scheduler")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# run a single simulation as part of the batch to fill a node
def run(sim):
    if args.task == 0:
        params["sim"] = sim
        params["beta"] = betas[sim]
        params["seed"] = random.randrange(1e9)
        file_name = "tutorial_run"+str(sim)+".txt"
        mc_lj(params, file_name=file_name)
        subprocess.call("~/feasst/build/bin/fst < " + file_name + " > tutorial_run"+str(sim)+".log", shell=True, executable='/bin/bash')
    else:
        subprocess.call("~/feasst/build/bin/rst < checkpoint" + str(sim) + ".txt", shell=True, executable='/bin/bash')

if __name__ == "__main__":
    if args.run_type == 0:
        mc_lj()
        subprocess.call("~/feasst/build/bin/fst < tutorial.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(params["num_sims"]) as pool:
            pool.starmap(run, zip(range(0, params["num_sims"])))
    elif args.run_type == 2:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt", shell=True, executable='/bin/bash')
    else:
        assert("unrecognized run_type")
