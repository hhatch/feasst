#!/usr/bin/env python

# # Description and Software Disclaimer

# Title: Example
#
# Software Description:
#
# Authors: Daniel W. Siderius, PhD; Harold W. Hatch, PhD
#
# ------------SOFTWARE DISCLAIMER AND REDISTRIBUTION CONDITIONS----------------
#
# This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to Title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This software is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgement if the software is used.
#
# This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
#
# ------------DECLARATION OF CONFLICTING INTERESTS----------------
#
# Certain commercially available items may be identified in this paper. This identification does not imply recommendation by NIST, nor does it imply that it is the best available for the purposes described.

import feasst
import pyfeasst

from FEASST_Henry_Coefficient_Rigid import *
from Uncertainty_Tools import *

#Relevant CODATA Constants
amu = 1.6605390400e-27 #kg/amu
kB = 1.380649030e-23 #J/K
Na = 6.022140760E23  #1/mol
h = 6.6260701500E-34  #J s

#Confidence Level
conf_level = 0.95

#Temperature
T = 300.

#Reference State Conditions
p = 1. # bar
V = kB * T / (p * 1.e5)  #molecular volume in m3

test_input = { "adsorptive": "data.C2",
               "adsorbent": "./LTA_replicate",
               "temperature": T,
               "rcut": 15.,
               "ncoeffs": 3,
               "trials": 1.e4,
               "pair_type": "LJCoulEwald",
               "tail_type": "LFS",
               "Ewald":{ "k2max": 27, "alpha": 6.0},
               "scale": {"active": True, "factor": 10000.}
             }


# Run the Henry-law Calculations
#Kcoeffs, Kcoeffs_var, beta = HenryCoefficient(debug=False,seed=835279,**test_input)
Kcoeffs, Kcoeffs_var, beta = Parallel_HenryCoefficient(nthreads=4,seed=835279,**test_input)

# Calculate deltaS from the Henry Law coefficients
deltaS_canonical = np.log(Kcoeffs[0]) - beta*Kcoeffs[1]/Kcoeffs[0]

# Estimate confidence intervals
CI_Kcoeffs, CI_deltaS_canonical = deltaS_Uncertainty_Canonical(Kcoeffs, Kcoeffs_var, beta, test_input, conf_level)

#print(test_input["adsorptive"])
#print('dedeltaS_canonical, CI_deltaS_canonical)

print(test_input["adsorptive"])
print(test_input["adsorbent"])
print('   deltaS_{ads}^{\infty}/R = '+str(deltaS_canonical))
print('     confidence interval   = '+str(CI_deltaS_canonical))

#NOTE: the "CI_deltaS_canonical" variable is the confidence interval, based on the chosen confidence level,
#  and assumption of linearized uncertainty propagation and the Welch-Satterthwaite approximation for the
#  number of degrees of freedom in the (derived) deltaS_canonical.

