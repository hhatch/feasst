*************************
README
*************************

.. .. image:: https://anaconda.org/hhatch/feasst/badges/installer/conda.svg
..     :target: https://conda.anaconda.org/hhatch

.. image:: https://travis-ci.com/hhatch/feasst.svg?branch=master
    :target: https://travis-ci.com/hhatch/feasst

.. image:: https://codecov.io/gh/hhatch/feasst/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/hhatch/feasst

.. image:: https://www.codefactor.io/repository/github/hhatch/feasst/badge/master
    :target: https://www.codefactor.io/repository/github/hhatch/feasst/overview/master
    :alt: CodeFactor

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, modular program to conduct molecular and particle-based
simulations with flat-histogram Monte Carlo methods.

.. note::

   Manuscript: https://doi.org/10.6028/jres.123.004

   Website: https://pages.nist.gov/feasst/

   Website DOI: https://doi.org/10.18434/M3S095

   Code repository: https://github.com/usnistgov/feasst

   Discussion list: https://groups.google.com/a/list.nist.gov/d/forum/feasst

Features
================================================================================

.. image:: dev/sphinx/feasst.png
   :target: https://pages.nist.gov/feasst
   :align: right

:doc:`/plugin/README` contains the list of features available.
Some features include but are not limited to the following:

Simulation techniques

* Wang-Landau Monte Carlo
* Transition-matrix Monte Carlo
* Metropolis Monte Carlo
* Mayer-sampling Monte Carlo
* Configurational bias

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* Grand canonical ensemble

Intermolecular interactions

* Hard spheres
* Lennard-Jones with Yukawa, LRC, force shift
* Patchy particles
* Charged interactions with the Ewald summation
* Cylindrical and slit pore confinement

Modern software

* Interface with C++ or as a Python module
* OpenMP parallelization
* Checkpointing to save and restart simulations

Usage
================================================================================

The following example Lennard-Jones Monte Carlo simulation may be found in `Tutorial <tutorial/tutorial.html>`_.

.. literalinclude:: tutorial/tutorial.py
   :language: py

Build from source code
=======================

Python install
----------------

* SWIG is required. Version 3.0.12 is recommended if your current SWIG version does not work properly.

* Anaconda with Python 3 is recommended.

* First, it is recommended to make a new python virtual environment.

.. code-block:: bash

    python3 -m venv ~/feasst_env
    source ~/feasst_env/bin/activate
    pip install jupyter matplotlib pandas scipy # for tutorials

* In the next step, CMake attempts to find the python libraries.
  But you may want to specify them manually.

.. code-block:: bash

    mkdir /path/to/feasst/build
    cd /path/to/feasst/build
    cmake -DUSE_SWIG=ON ..
    make _feasst -j$CPU_COUNT
    make install -j$CPU_COUNT
    python ../py/test.py  # optional test

* For manually setting the python path.

.. code-block:: bash

    cmake -DUSE_SWIG=ON -DSET_PYTHON_PATH=ON -DPYTHON_INCLUDE_DIR=/path/to/anaconda/include/python3.7m -DPYTHON_LIBRARY=/path/to/anaconda/lib/libpython3.7m.so ..

* You can deactivate and delete your python virtual environment as follows:

.. code-block:: bash

    deactivate
    rm -r ~/feasst_env  # if you'd like to reinstall or update

C++ install
----------------

First, install the C++ library.

.. code-block:: bash

    mkdir /path/to/feasst/build
    cd /path/to/feasst/build
    cmake ..
    make install -j$CPU_COUNT
    make test         # optional test

Then, compile the specific simulation you wish to run (e.g., tutorial).

.. code-block:: bash

    cd /path/to/feasst/tutorial/
    mkdir build; cd build
    cmake ..
    make
    ./tutorial

CMake defaults to install in the build directory.
But you can also specify the path as follows.

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

Later, when you build your tutorial executable, if your build directory is not `~/feasst/build`, then specify the path to the build directory as follows:

.. code-block:: bash

    cmake -DCMAKE_PREFIX_PATH=/path/to/install/dir ..

Troubleshooting install
------------------------

OS X
~~~~~~~~

* SWIG (from Homebrew) is likely version 4, which sometimes causes a SEGFAULT when trying to run feasst.
  Try SWIG version 3 instead.

* Sometimes CMake has trouble finding anaconda, and if you use SET_PYTHON_PATH described above, you may need to look out for the .dylib instead of .so

CentOS
~~~~~~~~~

CMake and SWIG versions are usually too old.
Try cmake3 instead of cmake.
Otherwise try SWIG 3.

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
