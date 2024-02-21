******************************
Python interface with pybind11
******************************

The pybind11 interface is in recent active development.

.. code-block:: bash

    sudo [brew/apt/yum/dnf] install pybind11
    cd /path/to/feasst/build
    cmake -DUSE_PYBIND11=ON ..
    make -j4
    export LD_LIBRARY_PATH="/path/to/feasst/build/:$LD_LIBRARY_PATH"
    pip install /path/to/feasst/

.. code-block:: py

    import feasst
    mc = feasst.MonteCarlo()
    text_input = """RandomMT19937 seed 123
    Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.fstprt"""
    for line in text_input.split('\n'):
        feasst.parse(mc, line)
