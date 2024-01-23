"""
Example of running feasst as a server and interacting with a python client.
"""

import subprocess
import multiprocessing
import time
import argparse
import socket
from pyfeasst import fstio

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--port', type=int, default=54321, help='server client interface port')
PARSER.add_argument('--buffer_size', type=int, default=1000, help='server client interface port')
PARSER.add_argument('--procs_per_node', type=int, default=1, help='number of processors')
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
PARAMS['prefix'] = 'listen'
PARAMS['script'] = __file__
PARAMS['sim_id_file'] = PARAMS['prefix']+ '_sim_ids.txt'
PARAMS['procs_per_sim'] = 1
PARAMS['num_sims'] = PARAMS['num_nodes']*PARAMS['procs_per_node']

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Listen port {port} buffer_size {buffer_size}
""".format(**params))

def server(params):
    script_file = params['prefix']+'_run.txt'
    write_feasst_script(params, script_file)
    subprocess.call(params['feasst_install']+'/bin/fst < '+script_file, shell=True, executable='/bin/bash')

def client(params):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(("localhost", params['port']))
    sock.send(b'Hello')
    #sock.send(b'EndListen')
    message = sock.recv(params['buffer_size'])
    print(str(message))
#    for x in range(0, 10000):
#        print("Step 1")
#        sock.send(b'Hello')
#        #sock.send(b'EndListen')
#        print("Step 2")
#        print(str(sock.recv(1000)))
#        print(x)
    

if __name__ == '__main__':
    proc1 = multiprocessing.Process(target=server, args=(PARAMS,))
    proc2 = multiprocessing.Process(target=client, args=(PARAMS,))
    proc1.start()
    time.sleep(0.1)
    proc2.start()
    proc1.join()
    proc2.join()
#    fstio.run_simulations(params=PARAMS,
#                          sim_node_dependent_params=None,
#                          write_feasst_script=write_feasst_script,
#                          post_process=None,
#                          queue_function=fstio.slurm_single_node,
#                          args=ARGS)
