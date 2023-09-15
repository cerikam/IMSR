# ******************************************************************************************************
#
#                                    EIRENE Write Burned Salt
#
#     Purpose: Uses SCALE Obiwan to read data from TRITON .f71 file and obtain burned salt atom
#              density vector for the previous depletion step.
#
# *******************************************************************************************************

#############################################
#   Write Burned Salt Atom Density Vector
#############################################

import sys, re
import os.path
import argparse
import subprocess
import salts_wf
import math

def get_burned_salt_atom_dens(f71_name):
    # Use SCALE obiwan to get data for specified fuel salt isotopes and in the specified order. Output is atom densities only (no string element names)
    output = subprocess.run(["/home/sigma//codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
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
            densities[nuclide] = float(data[2])  # The [1] here is what causes the code to return only the densities at position 1 of the f71 file
    be9 = densities['be-9']
    f19 = densities['f-19']
    li6 = densities['li-6']
    li7 = densities['li-7']
    u234 = densities['u-234']
    u235 = densities['u-235']
    u236 = densities['u-236']
    u238 = densities['u-238']

    burned_salt_dens = [li6, li7, be9, f19, u234, u235, u236, u238]   # Atom density vector for burned salt

    return burned_salt_dens