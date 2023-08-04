import sys
import subprocess
import math
import argparse
import random
import unittest

# define parameters of a pure component MC patchy trimer fluid
# see Fig 9 of http://dx.doi.org/10.1063/1.4918557
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/trimer.fstprt", "rwca": 2**(1./6.),
    "max_particles": 100, "min_particles": 0, "min_sweeps": 1e3, "mu": -1.375, "beta": 1/0.275,
    "trials_per": 1e5, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e5, "num_nodes": 1, "procs_per_node": 16, "script": __file__, "min_window_size": 5}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["trial_rigid_cluster_weight"] = 1./params['max_particles']

# write fst script to run a single simulation
def mc_trimer(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file trimer_lnpi.txt bounds_file trimer_bounds.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.25 min_size {min_window_size}
Checkpoint file_name trimer_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt} cutoff0_1 {rwca} cutoff1_1 {rwca}
NeighborCriteria energy_maximum -0.5 site_type0 0 site_type1 0
Potential EnergyMap EnergyMapNeighborCriteria neighbor_index 0 Model LennardJonesForceShift
RefPotential Model LennardJonesForceShift cutoff {rwca} VisitModel VisitModelCell min_length {rwca}
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25 particle_type 0
TrialRigidCluster weight {trial_rigid_cluster_weight} neighbor_index 0
Log trials_per_write {trials_per} file_name trimer[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-8

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min]
RemoveTrial name TrialAdd
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialTransfer weight 2 particle_type 0 reference_index 0 num_steps 8
RemoveAnalyze name Log
Log trials_per_write {trials_per} file_name trimer[sim_index].txt
Movie trials_per_write {trials_per} file_name trimer[sim_index].xyz
AnalyzeCluster trials_per_write {trials_per} file_name trimer_cluster[sim_index]_eq.txt multistate true stop_after_iteration 100
AnalyzeCluster trials_per_write {trials_per} file_name trimer_cluster[sim_index].txt multistate true start_after_iteration 100
Tune trials_per_write {trials_per} file_name trimer_tune[sim_index].txt multistate true stop_after_iteration 100
Energy trials_per_write {trials_per} file_name trimer_en[sim_index].txt multistate true start_after_iteration 100
CriteriaUpdater trials_per_update {trials_per}
CriteriaWriter trials_per_write {trials_per} file_name trimer_crit[sim_index].txt
""".format(**params))

# write slurm script to fill nodes with simulations
def slurm_queue():
    with open("slurm.txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N {num_nodes} -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 0 --task $SLURM_ARRAY_TASK_ID
if [ $? == 0 ]; then
  echo "Job is done"
  scancel $SLURM_ARRAY_JOB_ID
else
  echo "Job is terminating, to be restarted again"
fi
echo "Time is $(date)"
""".format(**params))

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_type', '-r', type=int, default=0, help="0: run, 1: submit to queue")
parser.add_argument('--task', type=int, default=0, help="input by slurm scheduler. If >0, restart from checkpoint.")
args = parser.parse_args()

# after the simulation is complete, perform some tests
class TestFlatHistogramTrimer(unittest.TestCase):
    def test(self):
        import numpy as np
        import pandas as pd

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "trimer_launch.txt"
        mc_trimer(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > trimer_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst trimer_checkpoint.fst", shell=True, executable='/bin/bash')
    if syscode == 0:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    elif args.run_type == 1:
        slurm_queue()
        subprocess.call("sbatch --array=0-10%1 slurm.txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    else:
        assert False  # unrecognized run_type
