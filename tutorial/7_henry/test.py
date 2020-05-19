#!/usr/bin/env python

# Title: Example Calculation of the Henry constant and infinite-dilution isosteric heat of
#   adsorption of TraPPE Ethane in LTA All-silica Zeolite
#
# Software Description: Calculation is accomplished via Widom Insertion to obtain N=1 state
#   spatial averages, using FEASST to generate the ghost insertions and perform energy calculations.
#   Calculation is based on the Henry coefficient calculation described in "Toward Rational Design of
#   Metal–Organic Frameworks for Sensing Applications: Efficient Calculation of Adsorption
#   Characteristics in Zero Loading Regime" by L. Sarkisov, <em>J Phys Chem C</em>, 116(4):3025-3033,
#   2012 (https://doi.org/10.1021/jp210633w)
#
# Requirements: This FEASST Script uses a specific branch and commit of the FEASST
#   source code, release tag v0.6.0 (https://github.com/usnistgov/feasst/releases/tag/v0.6.0),
#   as it utilizes tutorial-specific functions. Checkout the specific release of FEASST via:
#
#   git checkout v0.6.0
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
from scipy import stats

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
               "trials": 1.e6,
               "pair_type": "LJCoulEwald",
               "tail_type": "LFS",
               "Ewald":{ "k2max": 27, "alpha": 6.0},
               "scale": {"active": True, "factor": 10000.}
             }


# Run the Henry-law Calculations
#  Serial Version
#nthreads = 1
#Kcoeffs, Kcoeffs_var, beta = HenryCoefficient(debug=False,seed=835279,**test_input)
#  Parallel Version
nthreads = 4
Kcoeffs, Kcoeffs_var, beta = Parallel_HenryCoefficient(nthreads=nthreads,seed=835279,**test_input)

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
# print('mass:', cell_mass)
# print('vol: ', cell_volume)
rhoS = (cell_mass * amu) / (cell_volume * 1.e-30) #Skeletal Density in kg/m3
# print('density:', rhoS)

# Calculate the Henry Coefficient
Kh = Kcoeffs[0] * beta / rhoS
# Units:
#   beta = mol/kJ
#   rhoS = kg/m3
#   Kcoeffs[0] = dimensionless
#   Kh = (mol/kJ)/(kg/m3) = (mol/kg)(m3/1000J) = (mmol/kg)(1/Pa)
#     that is, millimoles adsorbate per kilogram adsorbent per Pascal
# Uncertainty in Kh
nu_Ki = int(test_input["trials"] * nthreads)  #degrees of freedom  in the MC integrals
a0 = (beta/rhoS)
var_Kh = (a0**2)*Kcoeffs_var[0]
coverage_factor = stats.t.ppf(1-(1.-conf_level)/2., nu_Ki-1)
CI_Kh = np.sqrt(var_Kh/float(nu_Ki-1)) * coverage_factor

# Calculate the infinite-dilution Isosteric Heat
Qst_inf = 1./beta + Kcoeffs[1]/Kcoeffs[0]  #FEASST native units = kJ/mol
#   Uncertainty in Qst_inf
a0 = Kcoeffs[1]/(Kcoeffs[0]**2)
a1 = 1./Kcoeffs[0]
var_Qst_inf = (a0**2)*Kcoeffs_var[0] + (a1**2)*Kcoeffs_var[1]

nu_eff = (var_Qst_inf**2)/(
    (a0**4)*(Kcoeffs_var[0]**2)/float(nu_Ki) + (a1**4)*(Kcoeffs_var[1]**2)/float(nu_Ki)
    )
nu_eff = int(nu_eff)
coverage_factor = stats.t.ppf(1-(1.-conf_level)/2., nu_eff-1)
CI_Qst_inf = np.sqrt(var_Qst_inf/float(nu_eff-1)) * coverage_factor

#NOTE: Uncertainty estimates are based on linear uncertainty propagation. Degrees of
# freedom are estimated from the sample count for direct measurements or using
# the Welch-Satterwaithe equation for derived quantities
# Uncertainty is assumed to derive exclusively from the statistical variation
# in the Monte Carlo averages. I.E., no error arises from the values of beta,
# rhoS, etc.

# Output to screen
print(test_input["adsorptive"])
print(test_input["adsorbent"])
print('   Henry Coefficient            = '+str(Kh)+' mmol/(kg Pa)')
print('     +/-                        = '+str(CI_Kh)+' mmol/(kg Pa)')
print('   Isosteric Heat (inf. dilute) = '+str(Qst_inf)+' kJ/mol')
print('     +/-                        = '+str(CI_Qst_inf)+' kJ/mol')
print(' +/- ranges are the 95 % confidence intervals estimated via linear uncertainty propagation')
