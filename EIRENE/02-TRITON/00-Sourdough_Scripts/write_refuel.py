# ******************************************************************************************************
#
#                                    EIRENE Write Refuel
#
#    Purpose: Creates atom density vector for refuel salt based on refuel enrichment and refuel amount.
#
# *******************************************************************************************************

import sys, re
import os.path
import argparse
import subprocess
import salts_wf
import math

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

def get_refuel_atom_dens(refuel_amount, refuel_enrichment):
    # Returns atom densities for refuel salt assuming UF4 molar percent of 5.0%
    # Assuming fuel matrix array of: [Li-6, Li-7, Be-9, F-19, U-234, U-235, U-236, U-238]
    #
    # Note: Function takes refuel_rate in cm^3/hour, and depletion step (dep_step) in days.
    #

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