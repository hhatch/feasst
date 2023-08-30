"""
Example Mayer-sampling simulation with a rigid coarse-grained mAb model
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyfeasst import fstio
from pyfeasst import physical_constants
from pyfeasst import coarse_grain_pdb
from pyfeasst import accumulator

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/plugin/chain/particle/cg7mab2.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--reference_sigma', type=float, default=4.5,
                    help='reference potential is a hard sphere of this size at the hinge')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e4),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibration')
PARSER.add_argument('--production_iterations', type=int, default=int(1e2),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=0.1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=0.1, help='hours until termination')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--prefix', type=str, default='cg7mab2_', help='prefix for all output file names')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--max_restarts', type=int, default=10, help='Number of restarts in queue')
PARSER.add_argument('--num_nodes', type=int, default=1, help='Number of nodes in queue')
PARSER.add_argument('--scratch', type=str, default=None,
                    help='Optionally write scheduled job to scratch/logname/jobid.')
PARSER.add_argument('--queue_flags', type=str, default="", help='extra flags for queue (e.g., for slurm, "-p queue")')
PARSER.add_argument('--node', type=int, default=0, help='node ID')
PARSER.add_argument('--queue_id', type=int, default=-1, help='If != -1, read args from file')
PARSER.add_argument('--queue_task', type=int, default=0, help='If > 0, restart from checkpoint')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['minutes'] = int(PARAMS['hours_terminate']*60) # minutes allocated on queue
PARAMS['hours_terminate'] = 0.99*PARAMS['hours_terminate'] - 0.0333 # terminate before queue
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length 500 particle_type0 /feasst/plugin/chain/particle/cg7mab2.fstprt \
    add_particles_of_type0 2 \
    group0 first first_particle_index 0
Potential Model HardSphere
RefPotential Model HardSphere sigma 0 sigma0 {reference_sigma} cutoff 0 cutoff0 {reference_sigma}
ThermoParams beta 1
MayerSampling num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
TrialTranslate new_only true reference_index 0 tunable_param 1 group first
TrialRotate new_only true reference_index 0 tunable_param 40
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}

# tune trial parameters
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}{sim}_b2_eq.txt
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.xyz
Tune
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name CriteriaWriter
RemoveAnalyze name Log
RemoveAnalyze name Movie

# production
CriteriaWriter trials_per_write {trials_per_iteration} file_name {prefix}{sim}_b2.txt
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}.xyz
MayerSampling num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    def b2(file_name):
        file1 = open(file_name, 'r')
        lines = file1.readlines()
        file1.close()
        exec('iprm=' + lines[0], globals())
        return iprm
    b2hs_ref = 2*np.pi*params['reference_sigma']**3/3 # reference HS in nm^3
    lpm_fac = 1e-24*physical_constants.AvogadroConstant().value()
    mlmg_fac = 1e-24*physical_constants.AvogadroConstant().value()/150000
    b2_overall = accumulator.Accumulator()
    for i in range(params['num_sims']):
        print('i', i)
        b2t = b2(params['prefix']+str(i)+'_b2.txt')
        print('b2', b2t['second_virial_ratio']*b2hs_ref, '+/-', b2t['second_virial_ratio_block_stdev']*b2hs_ref)
        print('b2', b2t['second_virial_ratio']*b2hs_ref*lpm_fac,
             '+/-', b2t['second_virial_ratio_block_stdev']*b2hs_ref*lpm_fac,
             'L/mol')
        b2_mlmg = b2t['second_virial_ratio']*b2hs_ref*mlmg_fac
        print('b2', b2_mlmg,
             '+/-', b2t['second_virial_ratio_block_stdev']*b2hs_ref*mlmg_fac,
             'mL/mg')
        b2_overall.add(b2_mlmg)
        assert np.abs(b2_mlmg - 0.0113) < 0.001
    print('b2 overall(ml/mg)', b2_overall.mean(), b2_overall.stdev()/np.sqrt(b2_overall.num_values()))

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          sim_node_dependent_params=None,
                          write_feasst_script=write_feasst_script,
                          post_process=post_process,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS)
