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

# Compute the cell mass and volume from the XYZ file
# NOTE: This is specific to SiO2 cells.
cell_mass = 0.
with open(test_input["adsorbent"]+'.xyz',mode='r') as f:
    lines = f.readlines()
    for i,line in enumerate(lines):
        if i == 1:
             cell = [float(x) for x in line.split()]
        elif i > 1:
            atom = line.split()[0]
            if atom == 'Si':
                cell_mass += 28.0855  #amu
            elif atom == 'O':
                cell_mass += 15.999   #amu
cell_volume = cell[0]*cell[1]*cell[2] #ang^3
print('mass:', cell_mass)
print('vol: ', cell_volume)
rhoS = (cell_mass * amu) / (cell_volume * 1.e-30) #Skeletal Density in kg/m3
print('density:', rhoS)

# Calculate the Henry Coefficient
Kh = Kcoeffs[0] * beta / rhoS
# Units:
#   beta = mol/kJ
#   rhoS = kg/m3
#   Kcoeffs[0] = dimensionless
#   Kh = (mol/kJ)/(kg/m3) = (mol/kg)(m3/1000J) = (mmol/kg)(1/Pa)
#     that is, millimoles adsorbate per kilogram adsorbent per Pascal
# Uncertainty in Kh


# Calculate the infinite-dilution Isosteric Heat
Qst_inf = 1./beta + Kcoeffs[1]/Kcoeffs[0]  #FEASST native units = kJ/mol
#   Uncertainty in Qst_inf


# Output to screen
print(test_input["adsorptive"])
print(test_input["adsorbent"])
print('   Henry Coefficient            = '+str(Kh)+' mmol/(kg Pa)')
#print('     confidence interval   = '+str(CI_deltaS_canonical))
print('   Isosteric Heat (inf. dilute) = '+str(Qst_inf)+' kJ/mol')
