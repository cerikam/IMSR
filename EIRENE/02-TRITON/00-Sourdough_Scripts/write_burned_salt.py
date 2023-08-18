# *************************************************************************************
#
#                      ----- EIRENE WRITE BURNED SALT -----
#
#     Purpose: This script uses the obiwan SCALE module to read data from a TRITON
#              .f71 file. The output is an atom density matrix for the burned fuel salt
#               from the previous depletion step in a TRITON simulation.
#
#     Note: This code will work for a standard .f71 TRITON file (no need to modify the
#           .f71 file to make it only one case)
#
# **************************************************************************************

import sys, re
import os.path
import argparse
import subprocess

def get_addnux4_nucs():

    addnux4 = ["h-1", "h-2", "h-3", "he-3", "he-4", "li-6", "li-7", "be-7", "be-9", "b-10", "b-11", "n-14 ", "n-15",
            "o-16", "o-17", "f-19", "na-23", "mg-24", "mg-25", "mg-26", "al-27", "si-28", "si-29", "si-30", "p-31",
            "s-32", "s-33", "s-34", "s-36", "cl-35", "cl-37", "ar-36", "ar-38", "ar-40", "k-39", "k-40", "k-41",
            "ca-40", "ca-42", "ca-43", "ca-44", "ca-46", "ca-48", "sc-45", "ti-46", "ti-47", "ti-48", "ti-49",
            "ti-50", "cr-50", "cr-52", "cr-53", "cr-54", "mn-55", "fe-54", "fe-56", "fe-57", "fe-58", "co-58",
            "co-59", "ni-58", "ni-59", "ni-60", "ni-61", "ni-62", "ni-64", "cu-63", "cu-65", "ga-69", "ga-71",
            "ge-70", "ge-72", "ge-73", "ge-74", "ge-76", "as-74", "as-75", "se-74", "se-76", "se-77", "se-78",
            "se-79", "se-80", "se-82", "br-79", "br-81", "kr-78", "kr-80", "kr-82", "kr-83", "kr-84", "kr-85",
            "kr-86", "rb-85", "rb-86", "rb-87", "sr-84", "sr-86", "sr-87", "sr-88", "sr-89", "sr-90", "y-89", 
            "y-90", "y-91", "zr-90", "zr-91", "zr-92", "zr-93", "zr-94", "zr-95", "zr-96", "nb-93", "nb-94", 
            "nb-95", "mo-92", "mo-94", "mo-95", "mo-96", "mo-97", "mo-98", "mo-99", "mo-100", "tc-99", "ru-96",
            "ru-98", "ru-99", "ru-100", "ru-101", "ru-102", "ru-103", "ru-104", "ru-105", "ru-106", "rh-103",
            "rh-105", "pd-102", "pd-104", "pd-105", "pd-106", "pd-107", "pd-108", "pd-110", "ag-107", "ag-109",
            "ag-111", "cd-106", "cd-108", "cd-110", "cd-111", "cd-112", "cd-113", "cd-114", "cd-116", "in-113",
            "in-115", "sn-112", "sn-113", "sn-114", "sn-115", "sn-116", "sn-117", "sn-118", "sn-119", "sn-120",
            "sn-122", "sn-123", "sn-124", "sn-125", "sn-126", "sb-121", "sb-123", "sb-124", "sb-125", "sb-126",
            "te-120", "te-122", "te-123", "te-124", "te-125", "te-126", "te-128", "te-130", "te-132", "i-127",
            "i-129", "i-130", "i-131", "i-135", "xe-123", "xe-124", "xe-126", "xe-128", "xe-129", "xe-130",
            "xe-131", "xe-132", "xe-133", "xe-134", "xe-135", "xe-136", "cs-133", "cs-134", "cs-135", "cs-136",
            "cs-137", "ba-130", "ba-132", "ba-133", "ba-134", "ba-135", "ba-136", "ba-137", "ba-138", "ba-140",
            "la-138", "la-139", "la-140", "ce-136", "ce-138", "ce-139", "ce-140", "ce-141", "ce-142", "ce-143",
            "ce-144", "pr-141", "pr-142", "pr-143", "nd-142", "nd-143", "nd-144", "nd-145", "nd-146", "nd-147",
            "nd-148", "nd-150", "pm-147", "pm-148", "pm-149", "pm-151", "sm-144", "sm-147", "sm-148", "sm-149",
            "sm-150", "sm-151", "sm-152", "sm-153", "sm-154", "eu-151", "eu-152", "eu-153", "eu-154", "eu-155",
            "eu-156", "eu-157", "gd-152", "gd-153", "gd-154", "gd-155", "gd-156", "gd-157", "gd-158", "gd-160",
            "tb-159", "tb-160", "dy-156", "dy-158", "dy-160", "dy-161", "dy-162", "dy-163", "dy-164", "ho-165",
            "er-162", "er-164", "er-166", "er-167", "er-168", "er-170", "lu-175", "lu-176", "hf-174", "hf-176",
            "hf-177", "hf-178", "hf-179", "hf-180", "ta-181", "ta-182", "w-182", "w-183", "w-184", "w-186", "re-185",
            "re-187", "ir-191", "ir-193", "au-197", "hg-196", "hg-198", "hg-199", "hg-200", "hg-201", "hg-202",
            "hg-204", "pb-204", "pb-206", "pb-207", "pb-208", "bi-209", "ra-223", "ra-224", "ra-225", "ra-226",
            "ac-225", "ac-226", "ac-227", "th-227", "th-228", "th-229", "th-230", "th-232", "th-233", "th-234",
            "pa-231", "pa-232", "pa-233", "u-232", "u-233", "u-234", "u-235", "u-236", "u-237", "u-238", "u-239",
            "u-240", "u-241", "np-235", "np-236", "np-237", "np-238", "np-239", "pu-236", "pu-237", "pu-238",
            "pu-239", "pu-240", "pu-241", "pu-242", "pu-243", "pu-244", "pu-246", "am-241", "am-242", "am-243",
            "am-244", "cm-241", "cm-242", "cm-243", "cm-244", "cm-245", "cm-246", "cm-247", "cm-248", "cm-249",
            "cm-250", "bk-249", "bk-250", "cf-249", "cf-250", "cf-251", "cf-252", "cf-253", "cf-254", "es-253",
            "es-254", "es-255", "co-58m", "ag-110m", "cd-115m", "te-127m", "te-129m", "pm-148m", "ho-166m",
            "am-242m", "am-244m"]

    return addnux4

def get_densities_from_f71(f71_name):
    
    # Use SCALE obiwan to get data from .f71 file
    output = subprocess.run(["./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
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
            densities[nuclide] = float(data[1])  # The [1] here is what causes the code to return only the densities at position 1 of the f71 file

    return densities

def get_burned_salt_atom_den(f71_name):
    
    # Extract atom density values for each isotope in fuel salt defined in list phermes_fuel

    phermes_fuel = ["be-9", "f-19", "li-6", "li-7", "u-234", "u-235", "u-236", "u-238"]
    min_dens = 1e-20
    densities = get_densities_from_f71(f71_name)

    # filter densities according to addnux4 and initial fuel composition
    nuclide_filter = list(set(get_addnux4_nucs()) | set(phermes_fuel))
    phermes_densities = list(filter(lambda d: d[0] in nuclide_filter, densities.items()))

    # Returns a list of atom density values in the order of the fuel materials listed in phermes_fuel
    atom_dens = []
    for nuclide, value in phermes_densities:
        for i in range(1):
            if value > min_dens:
                atom_dens.append(value)
    return atom_dens
