"""
Isothermal-isobaric ensemble Monte Carlo simulation of Lennard Jones particles.
"""

import argparse
import json
from pyfeasst import feasstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--fstprt', type=str, default='/feasst/particle/lj.fstprt',
                    help='FEASST particle definition')
PARSER.add_argument('--beta', type=float, default=1./0.9, help='inverse temperature')
PARSER.add_argument('--num_particles', type=int, default=500, help='number of particles')
PARSER.add_argument('--pressures', type=json.loads, default='{"pressure":[8.9429E-04, 2.6485E-03]}',
#NPT diverges from NVT when approaching the binodal.
#PARSER.add_argument('--pressures', type=json.loads, default='{"pressure":[8.9429E-04, 2.6485E-03, 4.3569E-03, 6.0193E-03, 7.6363E-03]}',
                    help='dictionary with a list of pressures to simulate')
PARSER.add_argument('--initial_cubic_side_length', type=int, default=20, help='cubic periodic boundary length')
PARSER.add_argument('--trials_per_iteration', type=int, default=int(1e5),
                    help='like cycles, but not necessary num_particles')
PARSER.add_argument('--equilibration_iterations', type=int, default=int(1e1),
                    help='number of iterations for equilibraiton')
PARSER.add_argument('--production_iterations', type=int, default=int(1e3),
                    help='number of iterations for production')
PARSER.add_argument('--hours_checkpoint', type=float, default=1, help='hours per checkpoint')
PARSER.add_argument('--hours_terminate', type=float, default=1, help='hours until termination')
PARSER.add_argument('--prefix', type=str, default='lj', help='prefix for all output file names')
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
PARAMS['procs_per_node'] = len(PARAMS['pressures']['pressure'])
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']
def sim_node_dependent_params(params):
    """ Set parameters that depent upon the sim or node here. """
    params['pressure'] = params['pressures']['pressure'][params['sim']]

def write_feasst_script(params, file_name):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(file_name, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
RandomMT19937 seed {seed}
Configuration cubic_side_length {initial_cubic_side_length} particle_type0 {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta {beta} chemical_potential -1
Metropolis
TrialTranslate tunable_param 2 tunable_target_acceptance 0.2
Checkpoint file_name {prefix}{sim}_checkpoint.fst num_hours {hours_checkpoint} num_hours_terminate {hours_terminate}
CheckEnergy trials_per_update {trials_per_iteration} tolerance 1e-4

# gcmc initialization
TrialAdd particle_type 0
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_init.txt
Tune
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
RemoveAnalyze name Log

# npt equilibration
ThermoParams beta {beta} pressure {pressure}
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {equilibration_iterations}
TrialVolume weight 0.1 tunable_param 0.2 tunable_target_acceptance 0.5
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}_eq.xyz
Density trials_per_write {trials_per_iteration} file_name {prefix}{sim}_density_eq.txt
Run until_criteria_complete true
RemoveModify name Tune
RemoveAnalyze name Log
RemoveAnalyze name Movie
RemoveAnalyze name Density

# npt production
Metropolis num_trials_per_iteration {trials_per_iteration} num_iterations_to_complete {production_iterations}
Log trials_per_write {trials_per_iteration} file_name {prefix}{sim}.txt
Movie trials_per_write {trials_per_iteration} file_name {prefix}{sim}.xyz start_after_iteration 1
Energy trials_per_write {trials_per_iteration} file_name {prefix}{sim}_en.txt
Density trials_per_write {trials_per_iteration} file_name {prefix}{sim}_density.txt
Volume trials_per_write {trials_per_iteration} file_name {prefix}{sim}_volume.txt
Run until_criteria_complete true
""".format(**params))

def post_process(params):
    """ Plot energy and compare with https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    ens = np.zeros(shape=(params['num_sims'], 2))
    rhos = np.zeros(shape=(params['num_sims'], 2))
    for sim in range(params['num_sims']):
        log = pd.read_csv(params['prefix']+str(sim)+'.txt')
        assert int(log['num_particles_of_type0'][0]) == params['num_particles']
        energy = pd.read_csv(params['prefix']+str(sim)+'_en.txt')
        ens[sim] = np.array([energy['average'][0],
                             energy['block_stdev'][0]])/params['num_particles']
        volume = pd.read_csv(params['prefix']+str(sim)+'_volume.txt', header=None)
        print('volume', volume)
        print('volume[0][0]', volume[0][0])
        print('volume[2][0]', volume[2][0])
        rhos[sim] = np.array([params['num_particles']/volume[0][0],
                              np.sqrt(params['num_particles'])*volume[2][0]/volume[0][0]])
        print('rhos[sim]', rhos[sim])
    # data from https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
    #rhos_srsw = [0.001, 0.003, 0.005, 0.007, 0.009]
    print('rhos', rhos)
    ens_srsw = [-9.9165E-03, -2.9787E-02]
    #ens_srsw = [-9.9165E-03, -2.9787E-02, -4.9771E-02, -6.9805E-02, -8.9936E-02]
    en_stds_srsw = [1.89E-05, 3.21E-05]
    plt.errorbar(params['pressures']['pressure'], ens_srsw, en_stds_srsw, fmt='+', label='SRSW')
    plt.errorbar(params['pressures']['pressure'], ens[:, 0], ens[:, 1], fmt='x', label='FEASST')
    plt.xlabel(r'$P$', fontsize=16)
    plt.ylabel(r'$U/N$', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(params['prefix']+'_energy.png', bbox_inches='tight', transparent='True')
    if len(ens_srsw) == params['num_sims']: # compare with srsw exactly
        for sim in range(params['num_sims']):
            diff = ens[sim][0] - ens_srsw[sim]
            print(diff)
            assert np.abs(diff) < 10*np.sqrt(ens[sim][1]**2 + en_stds_srsw[sim]**2)

if __name__ == '__main__':
    feasstio.run_simulations(params=PARAMS,
                             sim_node_dependent_params=sim_node_dependent_params,
                             write_feasst_script=write_feasst_script,
                             post_process=post_process,
                             queue_function=feasstio.slurm_single_node,
                             args=ARGS)
