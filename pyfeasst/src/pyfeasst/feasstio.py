"""
This module provides some utility input / output functions for use with the FEASST simulation program.
"""

import sys
import json
import subprocess
from multiprocessing import Pool
from itertools import repeat
import numpy as np

def vector3d_to_list(vec):
    """
    Converts a swig stl vector to python list

    >>> from pyfeasst import feasstio
    >>> feasstio.vector3d_to_list([[[[0]]]])
    [[[[0]]]]
    """
    lst = list()
    for _, vec1 in enumerate(vec):
        lst2 = list()
        for _, vec2 in enumerate(vec1):
            lst3 = list()
            for _, vec3 in enumerate(vec2):
                lst3.append(vec3)
            lst2.append(lst3)
        lst.append(lst2)
    return lst

def read_checkpoint(filename):
    """
    Return contents of checkpoint file as a string

    >>> from pyfeasst import feasstio
    >>> table = feasstio.read_checkpoint('../../tests/tutorial_0_table.txt')
    >>> table[:13]
    '6867 11 11 11'
    """
    with open (filename, "r") as myfile:
        checkpoint=myfile.readlines()
    assert(len(checkpoint) == 1)  # checkpoint files should have only one line
    return checkpoint[0]

def all_sims_complete(filename, num_sims):
    """
    Read filename and see if all sim ID's from [0, num_sims-1] are present (e.g., complete)

    >>> from pyfeasst import feasstio
    >>> all_sims_complete('../../tests/lj_sim_ids.txt', 8)
    True
    >>> all_sims_complete('../../tests/lj_sim_ids2.txt', 8)
    False
    """
    with open (filename, "r") as file1:
        lines = file1.read().splitlines()
    ids = map(str, list(range(num_sims)))
    for line in lines:
        ids = [i for i in ids if i not in line]
    if len(ids) == 0:
        return True
    return False

def run_simulations(params, run_function, post_process_function, queue_function, run_type, slurm_id, slurm_task):
    """
    Run a simulation either locally in the shell or queue on slurm

    :param dict params:
        Must have the following keys:
        sim_id_file: filename to write simulation id's for later checking of status,
        num_sims: number of simulations,
        prefix: prefix for all output file names,
        max_restarts: maximum number of restarts,
        num_nodes: number of nodes,
        node: node index.
    :param function run_function:
        The name of the function to run that represents a single simulation,
        and has the first argument that is the integer id of the simulation
        and the second arguement as the params.
    :param function post_process_function:
        The name of the function to post process all simulations once complete,
        and has the only argument as the params.
    :param function queue_function:
        The name of the function to queue one node and has the only argument as the params.
    :param int run_type:
        0: run on local shell, 1: submit to slurm, 2: post-process.
    :param int slurm_id:
        Input by slurm scheduler. If != -1, read args from file.
    :param int slurm_task:
        Input by slurm scheduler. If > 0, restart from checkpoint.
    """
    with open(params['sim_id_file'], 'w') as file1:
        file1.close() # clear file, then append sim id when complete
    if run_type == 0: # run directly
        if slurm_id != -1: # if run from SLURM
            if slurm_task == 0: # read param file if not checkpoint
                with open('lj_params'+str(slurm_id)+'.json', 'r') as file1:
                    params = json.load(file1)
        else:
            with open('lj_params.json', 'w') as file1:
                file1.write(json.dumps(params, indent=2))
        with Pool(params['num_sims']) as pool:
            codes = pool.starmap(run_function, zip(range(0, params['num_sims']), repeat(params)))
            if np.count_nonzero(codes) > 0:
                sys.exit(1)
    elif run_type == 1: # queue on SLURM
        slurm_id_file = params['prefix']+ '_slurm_ids.txt'
        with open(slurm_id_file, 'w') as file1:
            file1.close() # empty file contents
        for node in range(params['num_nodes']):
            params['node'] = node
            queue_function(params)
            subprocess.call("sbatch --array=0-" + str(params['max_restarts']) + "%1 " + params['prefix'] + "_slurm.txt | awk '{print $4}' >> " + slurm_id_file, shell=True, executable='/bin/bash')
            with open(slurm_id_file, 'r') as file1:
                slurm_id = file1.read().splitlines()[-1]
            with open('lj_params'+slurm_id+'.json', 'w') as file1:
                file1.write(json.dumps(params, indent=2))
    elif run_type == 2: # post process
        post_process_function(params)
    else:
        assert False  # unrecognized run_type

if __name__ == "__main__":
    import doctest
    doctest.testmod()
