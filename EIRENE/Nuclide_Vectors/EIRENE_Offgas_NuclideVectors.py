#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:47:29 2025

@author: sigma
"""

###############################################################################
#
#           Print Nuclide Atom Density Vectors as JSON File - EIRENE
#
#     For Noble Metals or Off-Gas system .f71 files!
#
###############################################################################

import json
import numpy as np
import os
import subprocess
import sys, re
from numpy import genfromtxt

SCALE_bin_path = '/home/sigma/codes/SCALE/SCALE-6.3.1/bin/'

f71_name = 'noble_metals.f71'     # Name of the .f71 file - change to desired name for noble gases/metals

# ***** File names: *****

# Off-gas system: noble_gases.f71

# Noble metal tracking: noble_metals.f71

enr = 3.5    # Refuel salt enrichment level

################################################################################

def read_fuel_salt_volume(iter):
    'Reads the fuel salt volume at the current depletion step.'
    
    # Change to directory for the specified refuel salt enrichment level:
    if enr == 3.5:
        dep_path = os.path.expanduser('~/EIRENE13/3.5/dep_step_{}'.format(iter))
    elif enr == 5:
        dep_path = os.path.expanduser('~/EIRENE13/5/dep_step_{}'.format(iter))
    elif enr == 6:
        dep_path = os.path.expanduser('~/EIRENE13/6/dep_step_{}'.format(iter))
    elif enr == 7:
        dep_path = os.path.expanduser('~/EIRENE13/7/dep_step_{}'.format(iter))
    elif enr == 9:
        dep_path = os.path.expanduser('~/EIRENE13/9/dep_step_{}'.format(iter))
    elif enr == 10:
        dep_path = os.path.expanduser('~/EIRENE13/10/dep_step_{}'.format(iter))
    elif enr == 15:
        dep_path = os.path.expanduser('~/EIRENE13/15/dep_step_{}'.format(iter))
    elif enr == 19.5:
        dep_path = os.path.expanduser('~/EIRENE13/19.5/dep_step_{}'.format(iter))
    elif enr == 19.75:
        dep_path = os.path.expanduser('~/EIRENE13/19.75/dep_step_{}'.format(iter))
    else:
        print('Invalid enrichment entry')
        
    os.chdir(dep_path)
    
    filename = 'Salt_volume_dep_step_{}.out'.format(iter)
    volume = genfromtxt(filename)
    
    return volume

def get_end_of_BOC_adens():
    '''Returns the atom density at the end of the completed BOC run.'''
    
    # Change to directory for the specified refuel salt enrichment level:
    if enr == 3.5:
        BOC_path = os.path.expanduser('~/EIRENE13/3.5/BOC')
    elif enr == 5:
        BOC_path = os.path.expanduser('~/EIRENE13/5/BOC')
    elif enr == 6:
        BOC_path = os.path.expanduser('~/EIRENE13/6/BOC')
    elif enr == 7:
        BOC_path = os.path.expanduser('~/EIRENE13/7/BOC')
    elif enr == 9:
        BOC_path = os.path.expanduser('~/EIRENE13/9/BOC')
    elif enr == 10:
        BOC_path = os.path.expanduser('~/EIRENE13/10/BOC')
    elif enr == 15:
        BOC_path = os.path.expanduser('~/EIRENE13/15/BOC')
    elif enr == 19.5:
        BOC_path = os.path.expanduser('~/EIRENE13/19.5/BOC')
    elif enr == 19.75:
        BOC_path = os.path.expanduser('~/EIRENE13/19.75/BOC')
    else:
        print('Invalid enrichment entry')
    
    os.chdir(BOC_path) 
    
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
    output = output.stdout.decode().split("\n")
    densities = {} # densities[nuclide] = (density at position 0 of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower() # convert to all lower cases
            num = int(dummy.group("num")) # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"   # for metastable isotopes
            else:
                nuclide = elem + "-" + str(num)
            densities[nuclide] = float(data[1])  # Position 0 is isotope names, 1 is data in noble metals/gases .f71 files
    
            BOC_sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
    
    return BOC_sorted_densities

def get_burned_salt_adens(iter):
    'Returns the atom density for the specified isotope in atoms/barn-cm.'
    # Change to directory for the specified refuel salt enrichment level:
    if enr == 3.5:
        dep_path = os.path.expanduser('~/EIRENE13/3.5/dep_step_{}'.format(iter))
    elif enr == 5:
        dep_path = os.path.expanduser('~/EIRENE13/5/dep_step_{}'.format(iter))
    elif enr == 6:
        dep_path = os.path.expanduser('~/EIRENE13/6/dep_step_{}'.format(iter))
    elif enr == 7:
        dep_path = os.path.expanduser('~/EIRENE13/7/dep_step_{}'.format(iter))
    elif enr == 9:
        dep_path = os.path.expanduser('~/EIRENE13/9/dep_step_{}'.format(iter))
    elif enr == 10:
        dep_path = os.path.expanduser('~/EIRENE13/10/dep_step_{}'.format(iter))
    elif enr == 15:
        dep_path = os.path.expanduser('~/EIRENE13/15/dep_step_{}'.format(iter))
    elif enr == 19.5:
        dep_path = os.path.expanduser('~/EIRENE13/19.5/dep_step_{}'.format(iter))
    elif enr == 19.75:
        dep_path = os.path.expanduser('~/EIRENE13/19.75/dep_step_{}'.format(iter))
    else:
        print('Invalid enrichment entry')
        
    os.chdir(dep_path)
    
    output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
    output = output.stdout.decode().split("\n")
    densities = {} # densities[nuclide] = (density at position 0 of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower() # convert to all lower cases
            num = int(dummy.group("num")) # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"   # for metastable isotopes
            else:
                nuclide = elem + "-" + str(num)
            densities[nuclide] = float(data[1])  # Position 0 is isotope names, 1 is BOC data, 2 is EOC data)
            
            sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
        
    return sorted_densities

file_index_forvol = np.arange(1, 146, 1)   # Create list of depletion step directories

file_index = np.arange(1, 147, 1)   # Add one to max value of file_index_forvol

fuel_salt_vols = []
for i in file_index_forvol:
    dep_vol = read_fuel_salt_volume(i)
    fuel_salt_vols = np.append(fuel_salt_vols, dep_vol)
    
d = {}

for i in range(len(file_index)):
    d[i] = {}
    if i == 0:
        d[i]['day'] = 0
        d[i]['saltvolume'] = 10680600.0
        d[i]['nuclide'] = 0    # Timezero adens for all nuclides is zero
    elif i == 1:
        d[i]['day'] = 7
        d[i]['saltvolume'] = 10680600.0
        d[i]['nuclide'] = get_end_of_BOC_adens()
    else:
        d[i]['day'] = i*7
        d[i]['saltvolume'] = fuel_salt_vols[i - 2]
        d[i]['nuclide'] = get_burned_salt_adens(i - 1)
            
########## Dump Data into JSON File: ##########

main_path = os.path.expanduser('~/EIRENE13/{}'.format(enr))
os.chdir(main_path)
        
with open("EIRENE_NobleMetal_Nuclides_{}.json".format(enr), "w") as outfile: 
    json.dump(d, outfile, indent=4)