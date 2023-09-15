# *****************************************************************************
#
#        EIRENE Find Volume of Previous Depletion Step and Use to Calculate
#        Height of Fuel Salt in Gas Plenum at Current Depletion Step Before
#        Adding in New Refuel
#
# ******************************************************************************

import sys, re
import os.path
import argparse
import subprocess
import salts_wf
import math


def add_refuel_volume(f71_name, refuel_amount, V0):
    # Function which returns h - the total fuel salt height in the gas plenum - after refuel has been added
    # Use SCALE obiwan to get burned salt volume data from .f71 file
    # V0 is the initial volume from the BOC core
    output = subprocess.run(["/home/sigma//codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
    output = output.stdout.decode().split("\n")
    for line in output:
        data = line.split(',')
        if "volume" in data[0]:
            vol = data[1:]
            volume = float(vol[1])
    total_vol = volume + refuel_amount    # Total volume of fuel salt (original V0 from prev. dep step + added refuel volume for current dep. step)
    if total_vol >= 2*V0:
        new_vol = total_vol/2                        # Half the total fuel salt volume
        new_refuel_vol = new_vol - V0                # In case there is any extra refuel that needs to be added to the gas plenum (for when total_vol > V0)
        h = (new_refuel_vol)/(math.pi*200**2)        # Height of refuel in gas plenum (cm)
    else:
        total_refuel_vol = total_vol - V0
        h = (total_refuel_vol)/(math.pi*200**2)
    return h
        
#test = add_refuel_volume("EIRENE_TRITON_1step.f71", 150000, 1551537.9)
#print(test)