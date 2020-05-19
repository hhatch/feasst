import numpy as np
#import copy
import os
import ase
import ase.io
import random
import string

#Physical Constants
amu = 1.6605390400e-27 #kg/amu
kB = 1.380649030e-23 #J/K
Na = 6.022140760E23  #1/mol
h = 6.6260701500E-34  #J s

def read_feasst_molecule(filename):
    with open(filename,mode='r') as f:
        data = f.read().splitlines()

    random_string = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(10)])
    tempfile = 'tempfile.'+random_string

    with open(tempfile,mode='w') as f:
        write_line = True
        Check_line = False
        lines_checked = 0
        for line in data:
            # Identify 'VMD Labels' area
            if line == 'VMD Labels':
                Check_line = True  # instruct wrapper to check subsequent lines
                write_line = False # turn off line writes
            # Check lines including and after 'VMD Labels'
            if Check_line:
                lines_checked += 1
                # Look for blank line after 'VMD Labels' area
                if lines_checked > 2 and line == '':  #minimum to check includes 'VMD Labels', blankline
                    Check_line = False #stop checking
                    lines_checked = 0 #reset counter
                    write_line = True #start writing with next line
            else:
                if write_line: f.write(line) ; f.write('\n')
                
    molecule = ase.io.read(tempfile, format='lammps-data')
    os.remove(tempfile)
    
    return molecule

def deBroglie(mass,T):
    Lambda = h / np.sqrt(2.e0 * np.pi * (mass*amu) * kB * T)
    return Lambda

def translational_S(T,molecule,V,N):
    mass = sum(molecule.get_masses())
    Lambda = deBroglie(mass,T)
    
    # Analytic
    S_trans = 1.5 + np.log(V/(Lambda**3)) -np.log(np.prod(range(1,N+1)))/float(N)
    # Stirling lnN! = NlnN - N
    #S_trans = 2.5 + np.log(V/(float(N)*Lambda**3))
    # Stirling lnN! = NlnN -N + 1
    #S_trans = 2.5 + np.log(V/(float(N)*Lambda**3)) - 1./float(N)
    # Stirling lnN! = (N+1/2)lnN - N + (1/2)ln(2Pi)
    #S_trans = 2.5 + np.log(V/(float(N)*Lambda**3)) -np.log(2.*np.pi*float(N))/(2.*float(N))
    
    return S_trans #divided by kB

def rotational_S(T,molecule,sigma):
    MoI = molecule.get_moments_of_inertia(vectors=False)
    #print(MoI)
    if MoI[0] < 1.e-20 and MoI[1] < 1.e-20 and MoI[2] < 1.e-20:
        # Case 1: Monatomic, zero moment of inertia
        S_rot = 0.
    elif MoI[0] < 1.e-20 and np.absolute(MoI[1]-MoI[2]) < 1.e-12:
        # Case 2: Two degenerate MoI
        I_conv = np.sqrt(MoI[1]*MoI[2])*amu*(1.e-10)**2
        #print(I_conv, 'converted')
        Qrot = 1./sigma
        Theta = (h**2)/(8. * (np.pi**2) * I_conv * kB)
        Qrot = Qrot * T/Theta
        S_rot = 1. + np.log(Qrot)
    else:
        # Case 3: Generalized Polyatomic
        Qrot = np.sqrt(np.pi)/sigma
        for I in MoI:
            I_conv = I*amu*(1.e-10)**2
            Theta = (h**2)/(8. * (np.pi**2) * I_conv * kB)
            Qrot = Qrot * np.sqrt(T/Theta)
        S_rot = 3./2. + np.log(Qrot)
    
    return S_rot  #divided by kB

def radius_of_gyration(input_mol):
    # Center of Mass
    COM = input_mol.get_center_of_mass()
    # Masses
    masses = input_mol.get_masses()
    # Positions
    r = input_mol.get_positions()
    # Radius of Gyration
    Rg = sum( [masses[i]*(np.linalg.norm(r[i]-COM)**2 ) for i in range(len(masses)) ] )
    Rg /= sum(masses)
    Rg = np.sqrt(Rg)
    
    return Rg
