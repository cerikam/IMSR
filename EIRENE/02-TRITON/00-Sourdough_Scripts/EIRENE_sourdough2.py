# *************************************************************************************************
#                             EIRENE Sourdough Script - Omnibus Edition
#
#   --- FIX THIS: Something is strange in that the values that this script spits out are slightly
#                 different than the original write_new_salt.py.
#
#    WORKING modified EIRENE_sourdough.py script which allows obiwan to run from a directory other
#    than the SCALE bin directory.
# *************************************************************************************************

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

#############################################
#   Write Refuel Salt Atom Density Vector
#############################################

# Establish total refuel amount added from chosen refuel rate

def get_total_refuel_amount(refuel_rate, dep_step):
    # Set depletion step size in days and use to calculate total refuel added from a given refuel amount.
    # Note: Entry of depletion step in days is required for this function.
    #       Entry of refuel rate in cubic centimeters/hour is required for this function.

    dep_step_hr = dep_step*24                    # Convert depletion step in days to hours
    refuel_amount = refuel_rate*dep_step_hr      # Total refuel amount added
    return refuel_amount

# *****************************
#  Get Refuel Weight Fractions
# *****************************

def saltmix(mU: float = 5) -> str:
    """Calculates salt mixture, assuming the MSRR salt is a melt of two salts,
     UF4 salt and a 66.6% LiF 33.3%BeF2 eutectic FLiBe.
    input: UF4 mol%
    output: salt name string"""

    mLi: float = (100.0 - mU) * 2.0 / 3.0
    mBe: float = mLi / 2.0
    mysalt = f'{mLi:5.3f}%LiF + {mBe:5.3f}%BeF2 + {mU:5.3f}%UF4'
    return mysalt

def get_refuel_wf(refuel_enrichment, UF4molpct):
    s = salts_wf.Salt(saltmix(UF4molpct), refuel_enrichment)
    iso_list = s.wf_gen()
    return iso_list

def get_refuel_den(tempK, refuel_enrichment, UF4molpct):
    s = salts_wf.Salt(saltmix(UF4molpct), refuel_enrichment)
    dens = float(s.densityK(tempK))
    return dens

# *************************************
# Convert Weight Frac. to Atom Dens.
# *************************************

def get_refuel_atom_dens(refuel_rate, dep_step, refuel_enrichment):
    # Returns atom densities for refuel salt assuming UF4 molar percent of 5.0%
    # Assuming fuel matrix array of: [Li-6, Li-7, Be-9, F-19, U-234, U-235, U-236, U-238]
    #
    # Note: Function takes refuel_rate in cm^3/hour, and depletion step (dep_step) in days.
    #

    refuel_amount = get_total_refuel_amount(refuel_rate, dep_step)

    refuel_wf = get_refuel_wf(refuel_enrichment, 5.0)
    refuel_dens = get_refuel_den(923.15, refuel_enrichment, 5.0)
    total_refuel_mass = refuel_amount*refuel_dens
    # Calculate mass of each isotope in the refuel:
    iso_mass = [i * total_refuel_mass for i in refuel_wf]
    # Calculate the number of atoms of each isotope in the refuel:
    atoms_Li6 = (iso_mass[0]*6.022e23)/6
    atoms_Li7 = (iso_mass[1]*6.022e23)/7
    atoms_Be9 = (iso_mass[2]*6.022e23)/9
    atoms_F19 = (iso_mass[3]*6.022e23)/19
    atoms_U234 = (iso_mass[4]*6.022e23)/234
    atoms_U235 = (iso_mass[5]*6.022e23)/235
    atoms_U236 = (iso_mass[6]*6.022e23)/236
    atoms_U238 = (iso_mass[7]*6.022e23)/238
    # List of atom values:
    iso_atoms = [atoms_Li6, atoms_Li7, atoms_Be9, atoms_F19, atoms_U234, atoms_U235, atoms_U236, atoms_U238]
    # Calculate atom density for each isotope
    atom_dens_cm3 = [i / refuel_amount for i in iso_atoms]   # Atom density of each isotope in atoms/cm3
    atom_dens_A = [i / 1e24 for i in atom_dens_cm3]          # Atom density of each isotope in atoms/cubic Angstrom
    return atom_dens_A

##################################################
#   Find height change from added refuel amount
##################################################

def add_refuel_to_core(refuel_rate, dep_step):
    # Finds the height change resulting from added refuel amount
    # Note: Function takes refuel_rate in cm^3/hour and depletion step (dep_step) in days.

    refuel_amount = get_total_refuel_amount(refuel_rate, dep_step)

    h = (refuel_amount)/(math.pi*200**2)   # Height of salt needed to add to gas plenum
    He_cyl_bottom = 440.5 + h              # New bottom z-dimension of He gas plenum
    return He_cyl_bottom

#######################################################################
#   Mix Refuel Salt with Burned Salt & Write New SCALE Material Input
#######################################################################

# ******************************
#  Extract burned salt volume
# ******************************

def get_burned_salt_volume(f71_name):
    # Use SCALE obiwan to get burned salt volume data from .f71 file
    output = subprocess.run(["/home/sigma//codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", f71_name], capture_output=True)
    output = output.stdout.decode().split("\n")
    for line in output:
        data = line.split(',')
        if "volume" in data[0]:
            vol = data[1:]
            volume = float(vol[1])
            return volume

# ***************************************************************************************************
#  Calculate number of atoms for each isotope in the burned fuel salt  (using atom density vector)
# ***************************************************************************************************

def burned_salt_atoms(f71_name):
    # Returns number of atoms for each isotope in burned fuel salt matrix
    # bs_volume = volume of the burned salt in cm^3
    bs_volume = get_burned_salt_volume(f71_name)     # Total burned salt volume in cm^3
    bs_adens = get_burned_salt_atom_dens(f71_name)   # Atom densities for each isotope in burned salt mat
    bs_vol_ang = bs_volume*1e24                      # Volume of burned salt in cubic Angstroms
    bs_atoms = [i * bs_vol_ang for i in bs_adens]    # Number of atoms for each isotope in burned salt mat
    return bs_atoms

# *********************************************************************************************
#  Calculate number of atoms for each isotope in the refuel salt (using atom density vector)
# *********************************************************************************************

def refuel_salt_atoms(refuel_rate, dep_step, refuel_enrichment):
    # Returns number of atoms for each isotope in refuel salt matrix
    # refuel_amount in cm^3
    # refuel_enrichment as a decimal (e.g. 0.03 for 3% U-235 enrichment)

    refuel_amount = get_total_refuel_amount(refuel_rate, dep_step)                  # Total refuel amount (cm^3) added for specified depletion step size

    refuel_adens = get_refuel_atom_dens(refuel_rate, dep_step, refuel_enrichment)   # Atom density for each isotope in refuel salt
    refuel_amount_ang = refuel_amount*1e24                                          # Refuel amount in cubic Angstroms                                
    refuel_atoms = [i * refuel_amount_ang for i in refuel_adens]
    return refuel_atoms

def mix_salt(f71_name, refuel_rate, dep_step, refuel_enrichment, bs_volume):
    # Mix burned salt and refuel salt

    refuel_amount = get_total_refuel_amount(refuel_rate, dep_step)                  # Total refuel amount (cm^3) added for specified depletion step size
    bs_volume = get_burned_salt_volume(f71_name)        # Volume of burned salt in cm^3

    bs_vol_ang = bs_volume*1e24                         # Volume of burned salt in cubic Angstroms
    refuel_amount_ang = refuel_amount*1e24              # Refuel amount in cubic Angstroms    
    tot_salt_vol = bs_vol_ang + refuel_amount_ang       # Total salt volume in cubic Angstroms
    bs_atoms = burned_salt_atoms(f71_name)   # Number of atoms for each isotope in burned salt mat

    # Splitting atom mat for each isotope in burned fuel:
    Li6_bs = bs_atoms[0]
    Li7_bs = bs_atoms[1]
    Be9_bs = bs_atoms[2]
    F19_bs = bs_atoms[3]
    U234_bs = bs_atoms[4]
    U235_bs = bs_atoms[5]
    U236_bs = bs_atoms[6]
    U238_bs = bs_atoms[7]

    refuel_atoms = refuel_salt_atoms(refuel_rate, dep_step, refuel_enrichment)    # Number of atoms for each isotope in refuel salt mat
    
    # Splitting refuel atom mat for each isotope in refuel:

    Li6_rf = refuel_atoms[0]
    Li7_rf = refuel_atoms[1]
    Be9_rf = refuel_atoms[2]
    F19_rf = refuel_atoms[3]
    U234_rf = refuel_atoms[4]
    U235_rf = refuel_atoms[5]
    U236_rf = refuel_atoms[6]
    U238_rf = refuel_atoms[7]

    # Finding total atoms for each isotope:

    Li6_tot = Li6_bs + Li6_rf
    Li7_tot = Li7_bs + Li7_rf
    Be9_tot = Be9_bs + Be9_rf
    F19_tot = F19_bs + F19_rf
    U234_tot = U234_bs + U234_rf
    U235_tot = U235_bs + U235_rf
    U236_tot = U236_bs + U236_rf
    U238_tot = U238_bs + U238_rf

    # Combining in one mat:
    mixed_atom_totals = [Li6_tot, Li7_tot, Be9_tot, F19_tot, U234_tot, U235_tot, U236_tot, U238_tot]

    # Finding atom densities of each from total salt volume (tot_salt_vol):
    mixed_atom_den = [i / tot_salt_vol for i in mixed_atom_totals]
    return mixed_atom_den

###############################################
#  Write new SCALE material composition block
###############################################

def write_scale_mat(f71_name, refuel_rate, dep_step, refuel_enrichment):
    # Returns a SCALE material composition block for the refuel salt mixed in with the burned salt

    bs_volume = get_burned_salt_volume(f71_name)     # Volume of burned salt in cm^3

    enrich_percent = refuel_enrichment*100    # Enrichment percent for refuel
    mixed_salt_aden = mix_salt(f71_name, refuel_rate, dep_step, refuel_enrichment, bs_volume)
    isotopes = ["li-6", "li-7", "be-9", "f-19", "u-234", "u-235", "u-236", "u-238"]
    scale_mat = "' Burned EIRENE fuel salt " + "and mixed refuel with enrichment of " + str(enrich_percent) + "%" + "\n"
    for i in range(8):
        scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], mixed_salt_aden[i])
    return scale_mat

#test = write_scale_mat("EIRENE03.f71",600,7,0.05)
#print(test)

#test = get_burned_salt_volume("EIRENE03.f71")
#print(test)
