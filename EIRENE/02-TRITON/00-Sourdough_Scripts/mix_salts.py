# ******************************************************************************************************
#
#                                    EIRENE Mix Salts
#
#    Purpose: Mixes burned salt and refuel salts and returns new SCALE material composition block.
#
# *******************************************************************************************************

import sys, re
import os.path
import argparse
import subprocess
import salts_wf
import math
import write_burned_salt as wb
import write_refuel as wr

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
    bs_adens = wb.get_burned_salt_atom_dens(f71_name)   # Atom densities for each isotope in burned salt mat
    bs_vol_ang = bs_volume*1e24                      # Volume of burned salt in cubic Angstroms
    bs_atoms = [i * bs_vol_ang for i in bs_adens]    # Number of atoms for each isotope in burned salt mat
    return bs_atoms

# *********************************************************************************************
#  Calculate number of atoms for each isotope in the refuel salt (using atom density vector)
# *********************************************************************************************

def refuel_salt_atoms(refuel_amount, refuel_enrichment):
    # Returns number of atoms for each isotope in refuel salt matrix
    # refuel_amount in cm^3
    # refuel_enrichment as a decimal (e.g. 0.03 for 3% U-235 enrichment)

    refuel_adens = wr.get_refuel_atom_dens(refuel_amount, refuel_enrichment)   # Atom density for each isotope in refuel salt
    refuel_amount_ang = refuel_amount*1e24                                          # Refuel amount in cubic Angstroms                                
    refuel_atoms = [i * refuel_amount_ang for i in refuel_adens]
    return refuel_atoms

def mix_salt(f71_name, refuel_amount, refuel_enrichment):
    # Mix burned salt and refuel salt

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

    refuel_atoms = refuel_salt_atoms(refuel_amount, refuel_enrichment)    # Number of atoms for each isotope in refuel salt mat
    
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

def write_scale_mat(f71_name, refuel_amount, refuel_enrichment):
    # Returns a SCALE material composition block for the refuel salt mixed in with the burned salt

    enrich_percent = refuel_enrichment*100    # Enrichment percent for refuel
    mixed_salt_aden = mix_salt(f71_name, refuel_amount, refuel_enrichment)
    isotopes = ["li-6", "li-7", "be-9", "f-19", "u-234", "u-235", "u-236", "u-238"]
    scale_mat = "' Burned EIRENE fuel salt " + "and mixed refuel with enrichment of " + str(enrich_percent) + "%" + "\n"
    for i in range(8):
        scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], mixed_salt_aden[i])
    return scale_mat