import sys
import subprocess
import argparse
import random
import unittest

# define parameters of a pure component NVT MC hard sphere simulation
params = {
    "cubic_box_length": 8, "fstprt": "/feasst/forcefield/atom.fstprt",
    "max_particles": 256, "min_particles": 0, "min_sweeps": 1e3, "mu": -2.352321,
    "trials_per": 1e6, "hours_per_adjust": 0.01, "hours_per_checkpoint": 1, "seed": random.randrange(int(1e9)), "num_hours": 5*24,
    "equilibration": 1e6, "num_nodes": 1, "procs_per_node": 32, "script": __file__, "min_window_size": 5}
params["num_minutes"] = round(params["num_hours"]*60)
params["hours_per_adjust"] = params["hours_per_adjust"]*params["procs_per_node"]
params["hours_per_checkpoint"] = params["hours_per_checkpoint"]*params["procs_per_node"]
params["num_hours_terminate"] = 0.95*params["num_hours"]*params["procs_per_node"]

# write fst script to run a single simulation
def mc_hs(params=params, file_name="launch.txt"):
    with open(file_name, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file hs_lnpi.txt bounds_file hs_bounds.txt num_adjust_per_write 10 min_window_size {min_window_size}
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.5 min_size {min_window_size}
Checkpoint file_name hs_checkpoint.fst num_hours {hours_per_checkpoint} num_hours_terminate {num_hours_terminate}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration cubic_box_length {cubic_box_length} particle_type0 {fstprt}
Potential Model HardSphere VisitModel VisitModelCell min_length 1
ThermoParams beta 1 chemical_potential {mu}
Metropolis
TrialTranslate weight 1 tunable_param 0.2 tunable_target_acceptance 0.25
Log trials_per_write {trials_per} file_name hs[sim_index].txt
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
Bias TransitionMatrix min_sweeps {min_sweeps} new_sweep 1
TrialTransfer weight 2 particle_type 0
Movie trials_per_write {trials_per} file_name hs[sim_index].xyz
Tune trials_per_write {trials_per} file_name hs_tune[sim_index].txt multistate true stop_after_iteration 100
#Energy trials_per_write {trials_per} file_name hs_en[sim_index].txt multistate true start_after_iteration 100
CriteriaUpdater trials_per_update {trials_per}
CriteriaWriter trials_per_write {trials_per} file_name hs_crit[sim_index].txt
PairDistribution trials_per_update 1000 trials_per_write {trials_per} \
  dr 0.025 file_name hs_gr[sim_index].txt multistate true multistate_aggregate false
#Run until_criteria_complete true
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
class TestFlatHistogramHS(unittest.TestCase):
    def test(self):
        # compare to EOS in SRSW: https://www.nist.gov/mml/csd/chemical-informatics-research-group/hard-sphere-thermodynamic-and-transport-properties
        import math
        import copy
        import numpy as np
        import pandas as pd
        from pyfeasst import macrostate_distribution

        lnpi = macrostate_distribution.MacrostateDistribution(file_name='hs_lnpi.txt')
        volume=8**3
        srsw=pd.read_csv('../test/data/stat_hs.csv')
        srsw=srsw[:6]
        pressure=list()
        lnpi_rw = copy.deepcopy(lnpi)
        for target_density in srsw['dens']:
            lnpi_rw.reweight_to_macrostate(target_macrostate=target_density*volume)
            pressure.append(-lnpi_rw.ln_prob().values[0]/volume)
        srsw['P_FST'] = pressure
        print(srsw[['dens', 'P_MC', 'P_FST', '+/-']])
        assert np.any(abs(srsw['P_MC'] - srsw['P_FST']) < 1e-3)

        # Use chemical potential from Carnahan-Starling to compare expected average density
        # http://www.sklogwiki.org/SklogWiki/index.php/Carnahan-Starling_equation_of_state
        rho = 0.1
        cubic_box_length=8
        eta = math.pi/6*rho
        betamu_ex = (8*eta-9*eta**2+3*eta**3)/(1-eta)**3
        betamu = betamu_ex + math.log(rho)
        lnpi_rw = lnpi.reweight(delta_beta_mu=betamu+2.352321)
        density=lnpi_rw.average_macrostate()/volume
        print('target_density', rho, 'density', density)
        assert abs(rho - density) < 5e-4

# run the simulation and, if complete, analyze.
def run():
    if args.task == 0:
        file_name = "hs_launch.txt"
        mc_hs(params, file_name=file_name)
        syscode = subprocess.call("../../../build/bin/fst < " + file_name + " > hs_launch.log", shell=True, executable='/bin/bash')
    else:
        syscode = subprocess.call("../../../build/bin/rst hs_checkpoint.fst", shell=True, executable='/bin/bash')
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
