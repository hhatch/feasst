"""
This module provides some utility input / output functions for use with the FEASST simulation program.
"""

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
