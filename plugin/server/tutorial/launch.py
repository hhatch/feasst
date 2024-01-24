"""
Example of running feasst as a server and interacting with a python client.
"""

import subprocess
import multiprocessing
import time
import argparse
import random
import socket
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--port', type=int, default=54321, help='server client interface port')
PARSER.add_argument('--buffer_size', type=int, default=1000, help='server client interface port')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
PARSER.add_argument('--seed', type=int, default=-1,
                    help='Random number generator seed. If -1, assign random seed to each sim.')
PARSER.add_argument('--run_type', '-r', type=int, default=0,
                    help='0: run, 1: submit to queue, 2: post-process')
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
PARAMS['prefix'] = 'listen'
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['procs_per_node']
if PARAMS['seed'] == -1:
    PARAMS['seed'] = random.randrange(int(1e9))

def client(params):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(("localhost", params['port']))
    for line in ['MonteCarlo',
                 'RandomMT19937 seed '+str(params['seed']),
                 'Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.fstprt',
                 'Potential Model LennardJones',
                 'ThermoParams beta 1 chemical_potential0 1',
                 'Metropolis',
                 'TrialTranslate',
                 'TrialAdd particle_type 0',
                 'Run until_num_particles 20',
                 'RemoveTrial name TrialAdd',
                 'Log output_file '+params['prefix']+'.csv',
                 'Run num_trials 10']:
        sock.send(bytes(line, 'utf-8'))
        message = sock.recv(params['buffer_size'])

if __name__ == '__main__':
    fstio.run_simulations(params=PARAMS,
                          queue_function=fstio.slurm_single_node,
                          args=ARGS,
                          client=client)
                          #clients=[client])
