#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:47:29 2025

@author: sigma
"""

###############################################################################
#
#           Print Nuclide Atom Density Vectors as JSON File - Th-EIRENE
#
###############################################################################

import json
import numpy as np
import os
import subprocess
import sys, re
from numpy import genfromtxt

SCALE_bin_path = '/home/sigma/codes/SCALE/SCALE-6.3.1/bin/'

f71_name = 'ThEIRENE.f71'     # Name of the .f71 file; change to ThEIRENE.f71 for ThEIRENE

################################################################################

def read_fuel_salt_volume(iter):
    'Reads the fuel salt volume at the current depletion step.'
    # Change to directory for the specified refuel salt enrichment level:
    dep_path = os.path.expanduser('~/ThEIRENE_Batch/dep_step_{}'.format(iter))
    os.chdir(dep_path)
    
    filename = 'Salt_volume_dep_step_{}.out'.format(iter)
    volume = genfromtxt(filename)
    
    return volume

def get_end_of_BOC_adens():
    '''Reads the new mixed fuel salt .f71 file produced by ORIGEN'''
    BOC_path = os.path.expanduser('~/ThEIRENE_Batch/BOC')        
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
            densities[nuclide] = float(data[2])  # Position 0 is isotope names, 1 is BOC data, 2 is EOC data
    
            BOC_sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
    
    return BOC_sorted_densities

def get_timezero_adens():
    'Returns the atom density at the BEGINNING of the BOC run for the specified isotope in atoms/barn-cm.'
    # Change to directory for the specified refuel salt enrichment level:
    BOC_path = os.path.expanduser('~/ThEIRENE_Batch/BOC')
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
            densities[nuclide] = float(data[1])  # Position 0 is isotope names, 1 is data)
    
            timezero_sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
    
    return timezero_sorted_densities

def get_burned_salt_adens(iter):
    'Returns the atom density for the specified isotope in atoms/barn-cm.'
    # Change to directory for the specified refuel salt enrichment level:
    dep_path = os.path.expanduser('~/ThEIRENE_Batch/dep_step_{}'.format(iter))
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
            densities[nuclide] = float(data[2])  # Position 0 is isotope names, 1 is BOC data, 2 is EOC data)
            
            sorted_densities = {k: v for k, v in sorted(densities.items(), key=lambda item: -item[1])}
        
    return sorted_densities

################################################################################

# ***** Testing Building a Dictionary ***** #

#file_index_forvol = np.arange(1, 366, 1)   # Create list of depletion step directories

#file_index = np.arange(1, 368, 1)

file_index_forvol = np.arange(1, 530, 1)   # Create list of depletion step directories

file_index = np.arange(1, 532, 1)

fuel_salt_vols = []
for i in file_index_forvol:
    dep_vol = read_fuel_salt_volume(i)
    fuel_salt_vols = np.append(fuel_salt_vols, dep_vol)

# Timezero adens:
    
adens_timezero = get_timezero_adens()

# End of BOC adens:
    
adens_BOCend = get_end_of_BOC_adens()

# Main depletion steps adens:

#for i in file_index:
#    adens_depmain = get_burned_salt_adens(i)

for i in range(len(file_index)):
    print(i)

d = {}

for i in range(len(file_index)):
    d[i] = {}
    if i == 0:
        d[i]['day'] = 0
        d[i]['saltvolume'] = 10680600.0
        d[i]['nuclide'] = get_timezero_adens()
    elif i == 1:
        d[i]['day'] = 7
        d[i]['saltvolume'] = 10680600.0
        d[i]['nuclide'] = get_end_of_BOC_adens()
    else:
        d[i]['day'] = i*7
        d[i]['saltvolume'] = fuel_salt_vols[i - 2]
        d[i]['nuclide'] = get_burned_salt_adens(i - 1)
        
########## Dump Data into JSON File: ##########

main_path = os.path.expanduser('~/ThEIRENE_Batch/')
os.chdir(main_path)
        
with open("10Year_ThEIRENE_FuelSalt_NuclideDensities.json", "w") as outfile: 
    json.dump(d, outfile, indent=4)