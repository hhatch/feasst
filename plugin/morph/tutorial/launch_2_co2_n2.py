import sys
import subprocess
import argparse
import random
import unittest
import pathlib
from pyfeasst import physical_constants

params = {
    "cubic_box_length": 30,
    "min_particles": 0,
    "fstprt0": "/feasst/forcefield/co2.fstprt",
    "fstprt1": "/feasst/forcefield/n2.fstprt",
    "temperature": 300, "max_particles": 10, "beta_mu": -15.24,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 1, "script": __file__}
params["ewald_alpha"] = 5.6/params["cubic_box_length"]
params["beta"] = 1./(params["temperature"]*physical_constants.MolarGasConstant(
).value()/1e3) # mol/kJ
params["mu"] = params["beta_mu"]/params["beta"]
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]
params["mu_init"] = -7
params["min_sweeps"]=1000
params["window_alpha"]=1.1
params["min_window_size"]=5

def mc_co2n2(params, file_name):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file co2n2_lnpin{node}.txt bounds_file co2n2_boundsn{node}.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha {window_alpha} min_size {min_window_size}
Checkpoint file_name co2n2_checkpointn{node}.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt0} particle_type1 {fstprt1}
Potential VisitModel Ewald alpha {ewald_alpha} kmax_squared 38
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential0 {mu_init} chemical_potential1 {mu_init}
Metropolis
TrialTranslate weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
Log trials_per_write {trials_per} file_name co2n2n{node}s[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1e-4

# gcmc initialization and nvt equilibration
TrialAdd particle_type 1
Run until_num_particles {max_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential0 {mu} chemical_potential1 {mu}
Metropolis
Run num_trials {equilibration}
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {max_particles} min {min_particles} soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 22 collect_flatness 20 min_collect_sweeps 18
TrialMorph particle_type0 0 particle_type_morph0 1
TrialMorph particle_type0 1 particle_type_morph0 0
RemoveAnalyze name Log
Log trials_per_write {trials_per} file_name co2n2n{node}s[sim_index].txt
Movie trials_per_write {trials_per} file_name co2n2n{node}s[sim_index].xyz
Tune trials_per_write {trials_per} file_name co2n2_tunen{node}s[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} file_name co2n2_enn{node}s[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per} file_name co2n2_critn{node}s[sim_index].txt
""".format(**params))

# write slurm script
def slurm_queue(node):
    params["node"]=node
    with open("slurm"+str(node)+".txt", "w") as myfile: myfile.write("""#!/bin/bash
#SBATCH -n {procs_per_node} -N 1 -t {num_minutes}:00 -o hostname_%j.out -e hostname_%j.out
echo "Running {script} ID $SLURM_JOB_ID on $(hostname) at $(date) in $PWD"
cd $PWD
export OMP_NUM_THREADS={procs_per_node}
python {script} --run_type 0 --task $SLURM_ARRAY_TASK_ID --node {node}
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
parser.add_argument('--node', type=int, default=0, help="break the job into multiple nodes.")
args = parser.parse_args()
params['node']=args.node

# after the simulation is complete, perform some tests
class TestFlatHistogramCO2N2(unittest.TestCase):
    def test(self):
        import numpy as np
        import pandas as pd
        lnpi=pd.read_csv('co2n2_lnpin0.txt')
        self.assertAlmostEqual(5.137717334901432, (np.exp(lnpi["ln_prob"]) * lnpi["state"]).sum(), delta=0.5)

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "co2n2_launch"+str(params["node"])+".txt"
        mc_co2n2(params=params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > co2n2_launch"+str(params['node'])+".log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst co2n2_checkpointn"+str(params['node'])+".fst", shell=True, executable='/bin/bash')
    if syscode == 0 and params['node'] == 1:
        unittest.main(argv=[''], verbosity=2, exit=False)
    return syscode

if __name__ == "__main__":
    if args.run_type == 0:
        syscode = run()
        if syscode != 0:
            sys.exit(1)
    elif args.run_type == 1:
        for node in range(params['num_nodes']):
            slurm_queue(node)
            subprocess.call("sbatch --array=0-2%1 slurm"+str(node)+".txt | awk '{print $4}' >> launch_ids.txt", shell=True, executable='/bin/bash')
    elif args.run_type == 2:
        unittest.main(argv=[''], verbosity=2, exit=False)
    else:
        assert False  # unrecognized run_type
