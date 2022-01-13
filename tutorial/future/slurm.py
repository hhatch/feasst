# queue a node with multiple simulations on an HPC with SLURM

import argparse
import subprocess
import random
from multiprocessing import Pool
import tutorial
import numpy as np

default_slurm_params = {
    "file_name": "slurm.txt",
    "num_nodes": 1,
    "procs_per_node": 4,
    "hours": 5*24}
default_slurm_params["num_sims"] = default_slurm_params["num_nodes"]*default_slurm_params["procs_per_node"]
betas = np.linspace(0.8, 1.2, num=default_slurm_params["num_sims"])

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
python slurm.py --run

echo "Time is $(date)"
""".format(**default_slurm_params))

def run(sim):
    params = tutorial.default_sim_params
    params["sim"] = sim
    params["beta"] = betas[sim]
    params["seed"] = random.randrange(1e9)
    file_name = "tutorial_run"+str(sim)+".txt"
    tutorial.mc_lj(params, file_name=file_name)
    subprocess.call("~/feasst/build/bin/fst < " + file_name + " > tutorial_run"+str(sim)+".log", shell=True, executable='/bin/bash')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', default=False, action='store_true')
    args = parser.parse_args()
    if args.run:
        with Pool(default_slurm_params["num_sims"]) as pool:
            pool.starmap(run, zip(range(0, default_slurm_params["num_sims"])))
    else:
        slurm_queue()
        subprocess.call("sbatch slurm.txt", shell=True, executable='/bin/bash')
