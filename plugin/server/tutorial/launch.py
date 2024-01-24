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


def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
Server port {port} buffer_size {buffer_size}
""".format(**params))

def server(params):
    script_file = params['prefix']+'_run.txt'
    write_feasst_script(params, script_file)
    subprocess.call(params['feasst_install']+'/bin/fst < '+script_file, shell=True, executable='/bin/bash')

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
        #print(str(message))

if __name__ == '__main__':
    proc1 = multiprocessing.Process(target=server, args=(PARAMS,))
    proc2 = multiprocessing.Process(target=client, args=(PARAMS,))
    proc1.start()
    time.sleep(0.1)
    proc2.start()
    proc1.join()
    proc2.join()
