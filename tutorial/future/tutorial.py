import subprocess
import numpy as np
import argparse
from multiprocessing import Pool
import random

parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: run single simulation on host, 1: run batch on host 2: submit match to scheduler")
args = parser.parse_args()

# define parameters a pure component NVT MC Lennard-Jones simulation
default_sim_params = {
    "seed": "time",
    "length": 8,
    "num_particles": 50,
    "fstprt": "/feasst/forcefield/lj.fstprt",
    "beta": 1.2,
    "steps_per": 1e5,
    "equilibration": 1e6,
    "production": 1e6,
    "sim": 0,
}

# write fst script to run a single simulation
def mc_lj(params=default_sim_params, file_name="tutorial.txt"):
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

# define scheduler params
default_slurm_params = {
    "file_name": "slurm.txt",
    "num_nodes": 1,
    "procs_per_node": 32,
    "hours": 5*24}
default_slurm_params["num_sims"] = default_slurm_params["num_nodes"]*default_slurm_params["procs_per_node"]

# set a simulation parameter to vary for for each processor
betas = np.linspace(0.8, 1.2, num=default_slurm_params["num_sims"])

# write slurm script
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node}
#SBATCH -N {num_nodes}
#SBATCH -t {hours}:00:00
#SBATCH -o hostname_%j.out
#SBATCH -e hostname_%j.out
echo "Running on host $(hostname)"
echo "Time is $(date)"
echo "Directory is $PWD"
echo "ID is $SLURM_JOB_ID"

cd $PWD
python tutorial.py --run_type 1

echo "Time is $(date)"
""".format(**default_slurm_params))

# run a single simulation as part of the batch to fill a node
def run(sim):
    params = default_sim_params
    params["sim"] = sim
    params["beta"] = betas[sim]
    params["seed"] = random.randrange(1e9)
    file_name = "tutorial_run"+str(sim)+".txt"
    mc_lj(params, file_name=file_name)
    subprocess.call("~/feasst/build/bin/fst < " + file_name + " > tutorial_run"+str(sim)+".log", shell=True, executable='/bin/bash')

if __name__ == "__main__":
    if args.run_type == 0:
        mc_lj()
        subprocess.call("~/feasst/build/bin/fst < tutorial.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 1:
        with Pool(default_slurm_params["num_sims"]) as pool:
            pool.starmap(run, zip(range(0, default_slurm_params["num_sims"])))
    elif args.run_type == 2:
        slurm_queue()
        subprocess.call("sbatch slurm.txt", shell=True, executable='/bin/bash')
    else:
        assert("unrecognized run_type")
