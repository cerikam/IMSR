#!/bin/env python3

import sys, re
import os.path
import argparse
import subprocess

###############################################################################

# This script is used to generate a SCALE composition block based on the
# nuclide densities stored in an f71 file.
#
# Nuclides are selected based on nuclides contained in the addnux4 set and
# based on a provided list of nuclides (e.g., the initial fuel composition).
# The f71 file should only have one case number on it. It is assumed the f71
# file was generated via blending using one of the other scripts in this
# directory. It cannot be applied to a standard TRITON-generated f71 file.

###############################################################################

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

    # get total mols via obiwan: initial actinide, final fission product
    output = subprocess.run(["{}/bin/obiwan".format(os.environ['SCALE']), "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
    output = output.stdout.decode().split("\n")

    densities = {} # densities[nuclide] = (density at position 0 of f71 file)
    skip = ["case", "step", "time", "power", "flux", "volume"]
    regexp = re.compile(r"(?P<elem>[a-zA-Z]+)(?P<num>\d+)(?P<meta>m)?")
    for line in output:
        data = line.split(',')
        if "case" in data[0]:
            case = set(data[1:])
            if len(case) != 1:
                sys.exit("Found more than one case on f71 file. Expected one case.")
        if data[0].strip() in skip:
            continue
        elif len(data) > 1:
            dummy = re.search(regexp, data[0].strip())
            elem = dummy.group("elem").lower() # convert to all lower cases
            num = int(dummy.group("num")) # to cut off leading zeros
            if dummy.group("meta"):
                nuclide = elem + "-" + str(num) + "m"
            else:
                nuclide = elem + "-" + str(num)
            densities[nuclide] = float(data[1])

    return densities


###############################################################################

def main():
    """
    Write core-average fuel composition
    """

    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument('-f71', '--f71-file', type=str,
                        dest='f71_name', required=True, metavar='',
                        help='core-average f71 file')

    parser.add_argument('-it', '--iteration', type=int,
                        default=0, dest='i_iter', metavar='',
                        help='Iteration number')

    args = parser.parse_args()

    phermes_fuel = ["u-234", "u-235", "u-236", "u-238", "c", "c-12", "c-13", "o", "o-16", "o-17"]
    mixture = 100
    min_dens = 1e-20
    i_iter = args.i_iter
    f71_name = os.path.abspath(args.f71_name)
    assert os.path.isfile(f71_name)

    if 'SCALE' not in os.environ:
        sys.exit("Export SCALE environmental variable")

    # get all densities from f71 file (first position on f71 file)
    print("Extracting densities from {}".format(f71_name))
    densities = get_densities_from_f71(f71_name)

    # filter densities according to addnux4 and initial fuel composition
    nuclide_filter = list(set(get_addnux4_nucs()) | set(phermes_fuel))
    phermes_densities = list(filter(lambda d: d[0] in nuclide_filter, densities.items()))

    # write to file in SCALE composition format
    out_filename = "core_average_density_{}.txt".format(i_iter)
    counter = 0
    with open(out_filename, "w") as FH:
        for nuclide, value in phermes_densities:
            if value > min_dens:
                counter += 1
                FH.write("{:10s} {:>3d} 0  {:>5e} 963.0 end \n".format(nuclide, mixture, value))
    print("Found {} nuclides on f71 file. Wrote {} into file {}.".format(len(densities), counter, out_filename))

###############################################################################
if __name__ == "__main__":
    main()
###############################################################################
