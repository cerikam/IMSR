# ******************************************************************************************
#
#                                  EIRENE Write New Salt
#       
#   Purpose: Mixes burned salt with refuel salt to create a new SCALE material input.
#            Returns atom densitiy vector for resulting mixed salt.
#
# *******************************************************************************************

import write_burned_salt as ws
import write_refuel as wr

def burned_salt_atoms(f71_name, bs_volume=1551567.8):
    # Returns number of atoms for each isotope in burned fuel salt matrix
    # bs_volume = volume of the burned salt in cm^3
    bs_adens = ws.get_burned_salt_atom_den(f71_name)   # Atom densities for each isotope in burned salt mat
    bs_vol_ang = bs_volume*1e24                         # Volume of burned salt in cubic Angstroms
    bs_atoms = [i * bs_vol_ang for i in bs_adens]    # Number of atoms for each isotope in burned salt mat
    return bs_atoms

def refuel_salt_atoms(refuel_amount, refuel_enrichment):
    # Returns number of atoms for each isotope in refuel salt matrix
    # refuel_amount in cm^3
    # refuel_enrichment as a decimal (e.g. 0.03 for 3% U-235 enrichment)
    refuel_adens = wr.get_refuel_atom_dens(refuel_amount, refuel_enrichment)   # Atom density for each isotope in refuel salt
    refuel_amount_ang = refuel_amount*1e24                                     # Refuel amount in cubic Angstroms                                
    refuel_atoms = [i * refuel_amount_ang for i in refuel_adens]
    return refuel_atoms

def mix_salt(f71_name, refuel_amount, refuel_enrichment, bs_volume=1551567.8):
    # Mix burned salt and refuel salt
    bs_vol_ang = bs_volume*1e24                         # Volume of burned salt in cubic Angstroms
    refuel_amount_ang = refuel_amount*1e24              # Refuel amount in cubic Angstroms    
    tot_salt_vol = bs_vol_ang + refuel_amount_ang       # Total salt volume in cubic Angstroms
    bs_atoms = burned_salt_atoms(f71_name, 1551567.8)   # Number of atoms for each isotope in burned salt mat
    # Splitting atom mat for each isotope in burned fuel:
    Be9_bs = bs_atoms[0]
    F19_bs = bs_atoms[1]
    Li6_bs = bs_atoms[2]
    Li7_bs = bs_atoms[3]
    U234_bs = bs_atoms[4]
    U235_bs = bs_atoms[5]
    U236_bs = bs_atoms[6]
    U238_bs = bs_atoms[7]
    #  Note: bs_atoms is in the order of:  phermes_fuel = ["be-9", "f-19", "li-6", "li-7", "u-234", "u-235", "u-236", "u-238"]
    refuel_atoms = refuel_salt_atoms(refuel_amount, refuel_enrichment)    # Number of atoms for each isotope in refuel salt mat
    # Splitting refuel atom mat for each isotope in refuel:
    Be9_rf = refuel_atoms[2]
    F19_rf = refuel_atoms[3]
    Li6_rf = refuel_atoms[0]
    Li7_rf = refuel_atoms[1]
    U234_rf = refuel_atoms[4]
    U235_rf = refuel_atoms[5]
    U236_rf = refuel_atoms[6]
    U238_rf = refuel_atoms[7]
    #  Note: refuel_atoms is in the order of: Li-6, Li-7, Be-9, F-19, U-234, U-235, U-236, U-238
    # Finding total atoms for each isotope:
    Be9_tot = Be9_bs + Be9_rf
    F19_tot = F19_bs + F19_rf
    Li6_tot = Li6_bs + Li6_rf
    Li7_tot = Li7_bs + Li7_rf
    U234_tot = U234_bs + U234_rf
    U235_tot = U235_bs + U235_rf
    U236_tot = U236_bs + U236_rf
    U238_tot = U238_bs + U238_rf
    # Combining in one mat:
    mixed_atom_totals = [Be9_tot, F19_tot, Li6_tot, Li7_tot, U234_tot, U235_tot, U236_tot, U238_tot]
    # Finding atom densities of each from total salt volume (tot_salt_vol):
    mixed_atom_den = [i / tot_salt_vol for i in mixed_atom_totals]
    return mixed_atom_den

def write_scale_mat(f71_name, refuel_amount, refuel_enrichment, bs_volume=1551567.8):
    # Returns a SCALE material composition block for the refuel salt mixed in with the burned salt
    enrich_percent = refuel_enrichment*100    # Enrichment percent for refuel
    mixed_salt_aden = mix_salt(f71_name, refuel_amount, refuel_enrichment, 1551567.8)
    isotopes = ["be-9", "f-19", "li-6", "li-7", "u-234", "u-235", "u-236", "u-238"]
    scale_mat = "' Burned EIRENE fuel salt " + "and mixed refuel with enrichment of " + str(enrich_percent) + "%" + "\n"
    for i in range(8):
        scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], mixed_salt_aden[i])
    return scale_mat
