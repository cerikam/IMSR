# ***************************************************************************************************
#
#                        EIRENE Sourdough Submit and Analyze KENO Runs
#
#   By: C. Erika Moss and Dr. Ondrej Chvala
#
# ****************************************************************************************************

import time
import os
import shutil
import salts_wf
import initialize_BOC
from initialize_BOC import BOC_core
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import math
import sys, re
import subprocess

#####################################################
#     Write KENO Decks Based on Refuel Parameters
#####################################################

class Refuel_KENO_deck(object):
    'Create parallel KENO decks for EIRENE for different refuel volume additions.'
    def __init__(self):
        # ************************
        #   Set Material Params.
        # ************************
        self.fs_tempK:float = 923.15     # Fuel salt temperature in Kelvin
        self.mat_tempK:float = 923.15    # Material temperature in Kelvin
        self.gr_tempK:float = 923.15     # Graphite temperature in Kelvin
        self.gr_dens:float = 2.3         # Graphite density in g/cc
        self.ss316_dens:float = 8.03     # SS316 steel density in g/cc
        self.hast_dens:float = 8.89      # Hastelloy N density in g/cc
        self.he_dens:float = 0.0001785   # Helium density in g/cc
        self.UF4molpct:float = 5.0       # UF4 mol % in FLiBe-U
        self.Uenrichpct:float = 2.650    # Initial Uranium enrichment %
        self.Uenrichfrac:float = 0.0265  # Initial Uranium enrichment fraction

        # ************************
        #    Set Refuel Params.
        # ************************
        self.renrich = 0.03                                # Refuel enrichment fraction
        self.rvols = [5000, 8000, 10000, 12000, 14000]     # Refuel volumes in cm^3
        self.V0 = 1551537.9                                # Total fuel salt volume for BOC core in cm^3

        # ************************
        #    Set KENO Params.
        # ************************
        self.npg:float = 40000           # Number per generation (npg)
        self.nsk:float = 50              # Number of skipped generations (nsk)
        self.gen:float = 5000            # Number of generations (gen)
        self.sig:float = 50e-5           # (sig)
        
        # ************************
        #    Running Parameters
        # ************************
        self.queue:str = 'fill'                     # NECluster queue
        self.ppn:int = 64                           # ppn core count
        self.deck_name:str = 'EIRENE.inp'           # KENO input file name
        self.f71_name:str = 'EIRENE.f71'            # Name of .f71 TRITON file to read
        self.qsub_name:str = 'runEIRENE-Scale.sh'   # Name of each qsub file
        self.qsub_path:str = os.path.expanduser('~/EIRENE/Scripts/Sourdough/runEIRENE-Scale.sh')  # Full path to the qsub script

        # *************************
        #    Data Output Params.
        # *************************

        self.crit_refuel:float = None               # Placeholder for critical refuel amount in cm^3

        ########################################################
        #  Write the Burned Salt from Previous Depletion Step
        ########################################################

    def get_burned_salt_atom_dens(self):
        # Use SCALE obiwan to get data for specified fuel salt isotopes and in the specified order. Output is atom densities only (no string element names)
        output = subprocess.run(["/home/sigma/codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", self.f71_name], capture_output=True)
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

    ################################
    #     Write the Refuel Salt
    ################################

    def saltmix(self, mU: float = 5) -> str:
        """Calculates salt mixture, assuming the MSRR salt is a melt of two salts,
        UF4 salt and a 66.6% LiF 33.3%BeF2 eutectic FLiBe.
        input: UF4 mol%
        output: salt name string"""

        mLi: float = (100.0 - mU) * 2.0 / 3.0
        mBe: float = mLi / 2.0
        mysalt = f'{mLi:5.3f}%LiF + {mBe:5.3f}%BeF2 + {mU:5.3f}%UF4'
        if True:
            print("Salt: ", mysalt)
            print("Molar percent sum: ", mBe + mLi + mU)
        return mysalt
    
    def get_refuel_wf(self):
        'Returns the weight fractions of each isotope component of the refuel salt.'
        s = salts_wf.Salt(self.saltmix(self.UF4molpct), self.renrich)
        iso_list = s.wf_gen()
        return iso_list

    def get_refuel_den(self):
        'Returns the density of the refuel salt.'
        s = salts_wf.Salt(self.saltmix(self.UF4molpct), self.renrich)
        dens = float(s.densityK(self.fs_tempK))
        return dens
    
    def get_refuel_atom_dens(self, rvol):
        'Returns the atom density vector of the refuel salt based on the weight fractions of each isotope.'
        refuel_wf = self.get_refuel_wf()
        refuel_dens = self.get_refuel_den()
        total_refuel_mass = rvol*refuel_dens
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
        atom_dens_cm3 = [i / rvol for i in iso_atoms]   # Atom density of each isotope in atoms/cm3
        atom_dens_A = [i / 1e24 for i in atom_dens_cm3]          # Atom density of each isotope in atoms/cubic Angstrom
        return atom_dens_A

        ######################################################################
        #  Mix Refuel in with Burned Salt and Write New SCALE Material Input
        ######################################################################

    def get_burned_salt_volume(self):
        'Returns the total volume of salt in the burned core from previous depletion step.'
        output = subprocess.run(["/home/sigma/codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", self.f71_name], capture_output=True)
        output = output.stdout.decode().split("\n")
        for line in output:
            data = line.split(',')
            if "volume" in data[0]:
                vol = data[1:]
                volume = float(vol[1])
                return volume

    def burned_salt_atoms(self):
        'Returns the total number of atoms for each isotope in the burned fuel salt matrix from the previous depletion step.'
        bs_volume = self.get_burned_salt_volume()     # Total burned salt volume in cm^3
        bs_adens = self.get_burned_salt_atom_dens()   # Atom densities for each isotope in burned salt mat
        bs_vol_ang = bs_volume*1e24                      # Volume of burned salt in cubic Angstroms
        bs_atoms = [i * bs_vol_ang for i in bs_adens]    # Number of atoms for each isotope in burned salt mat
        return bs_atoms
    

    def refuel_salt_atoms(self, rvol):
        'Returns the total number of atoms for each isotope in the refuel salt matrix.'
        refuel_adens = self.get_refuel_atom_dens(rvol)   # Atom density for each isotope in refuel salt
        refuel_amount_ang = rvol*1e24                                          # Refuel amount in cubic Angstroms                                
        refuel_atoms = [i * refuel_amount_ang for i in refuel_adens]
        return refuel_atoms
    
    def mix_salt(self, rvol):
        'Mix the burned salt from the previous step with the refuel salt.'
        bs_volume = self.get_burned_salt_volume()        # Volume of burned salt in cm^3

        bs_vol_ang = bs_volume*1e24                         # Volume of burned salt in cubic Angstroms
        refuel_amount_ang = rvol*1e24              # Refuel amount in cubic Angstroms    
        tot_salt_vol = bs_vol_ang + refuel_amount_ang       # Total salt volume in cubic Angstroms
        bs_atoms = self.burned_salt_atoms()   # Number of atoms for each isotope in burned salt mat

        # Splitting atom mat for each isotope in burned fuel:
        Li6_bs = bs_atoms[0]
        Li7_bs = bs_atoms[1]
        Be9_bs = bs_atoms[2]
        F19_bs = bs_atoms[3]
        U234_bs = bs_atoms[4]
        U235_bs = bs_atoms[5]
        U236_bs = bs_atoms[6]
        U238_bs = bs_atoms[7]

        refuel_atoms = self.refuel_salt_atoms(rvol)    # Number of atoms for each isotope in refuel salt mat
        
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

    def write_scale_mat(self, rvol):
        'Returns a SCALE material composition block for the refuel salt mixed in with the burned salt.'
        enrich_percent = self.renrich*100    # Enrichment percent for refuel
        mixed_salt_aden = self.mix_salt(rvol)
        isotopes = ["li-6", "li-7", "be-9", "f-19", "u-234", "u-235", "u-236", "u-238"]
        scale_mat = "' Burned EIRENE fuel salt " + "and mixed refuel with enrichment of " + str(enrich_percent) + "%" + "\n"
        for i in range(8):
            scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], mixed_salt_aden[i])
        return scale_mat
    
    #################################
    #        Add in Refuel
    #################################

    def add_refuel_volume(self, rvol):
        'Returns h, the total fuel salt height in the gas plenum after refuel salt has been added.'
        'The V0 variable is the initial volume from the BOC core.'
        output = subprocess.run(["/home/sigma/codes/SCALE/SCALE-6.3.1/bin/./obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", self.f71_name], capture_output=True)
        output = output.stdout.decode().split("\n")
        for line in output:
            data = line.split(',')
            if "volume" in data[0]:
                vol = data[1:]
                volume = float(vol[1])
        total_vol = volume + rvol   # Total volume of fuel salt (original V0 from prev. dep step + added refuel volume for current dep. step)
        # Check if the volume has doubled:
        if total_vol >= 2*self.V0:
            new_vol = total_vol/2                        # Half the total fuel salt volume
            new_refuel_vol = new_vol - self.V0                # In case there is any extra refuel that needs to be added to the gas plenum (for when total_vol > V0)
            h = (new_refuel_vol)/(math.pi*200**2)        # Height of refuel in gas plenum (cm)
        else:
            total_refuel_vol = total_vol - self.V0
            h = (total_refuel_vol)/(math.pi*200**2)
        return h
    
    ########################################
    #      Write Parallel KENO Decks
    ########################################

    def write_KENO_decks(self):
        'Creates a directory for each refuel amount and writes a KENO deck there.'
        for x in range(len(self.rvols)):
            scale_fuel = self.write_scale_mat(self.rvols[x])
            h = 440.5 + self.add_refuel_volume(self.rvols[x])
            deckpath = f'refuel_{self.rvols[x]:5.01f}'
            if not os.path.isdir(deckpath):
                os.mkdir(deckpath)
            os.chdir(deckpath)

            keno_deck = f'''=csas6 parm=(   )
EIRENE SCALE/CSAS model, UF4 mol% = {self.UF4molpct} and refuel enrichment {self.renrich}
ce_v7.1

read comp
' FLIBe-U fuel salt
{scale_fuel}

' Hastelloy N
wtptHastelloy 2 8.89 5
         28000 71.0
         24000 7.0
         26000 5.0
         14000 1.0
         42000 16.0
         1.0 923.15
         28058 67.6 28060 32.4
         24052 100.0
         26056 100.0
         14028 100.0
         42092 14.65 42094 9.19 42095 15.87 42096 16.67 42097 9.58 42098 24.29 42100 9.75 end
' SS Shutdown Rods
    ss316 3 den=2.7 1.0 923.15 end
' Graphite
   graphite 4 den=1.84 1.0 923.15 end
' Stainless Steel SS316
   ss316 5 den=8.030000 1.0 923.15 end
' Helium gas
   he 6 den=0.0001785 1.0 923.15 end

end comp

read parameters
 npg=5000 nsk=50 gen=230 sig=10e-5
 htm=no
 fdn=no
 pms=no
 pmv=no
end parameters

read geometry
' ****************************************************************************
'                             CONTROL RODS
' ****************************************************************************

' Verify control rod height

' Unit 1: Control rod unit cell (cylinder 1 = shutdown rod; cylinder 2 = guide tube)
unit 1
' Shutdown rod (Gd)
  cylinder 1 7.0 415.0 130
' Guide tube
  cylinder 2 7.2 415.0 130
' Fuel region beneath
  cylinder 3 7.2 130 0.0
' Graphite hexprism
  hexprism 4 13.076 415.0 0.0
  media 3 1 1
  media 2 1 -1 2
  media 1 1 3
  media 4 1 4 -1 -2 -3
  boundary 4

' ****************************************************************************
'                        REGION 1 FUEL ASSEMBLIES
' ****************************************************************************

' ----- Small molten salt fuel channel -----

' Unit 2: Region 1 Fuel Assembly - Single Cylindrical Unit
unit 2
' Single small fuel channel
  cylinder 1 1.157 415.0 0.0
' Graphite hexprism to contain fuel channel
  hexprism 2 2.6 415.0 0.0
  media 1 1 1
  media 4 1 -1 2
  boundary 2

' Unit 3: Region 1 Fuel Assembly - Bare Graphite Hexprism
unit 3
  hexprism 1 2.6 415.0 0.0
  media 4 1 1
  boundary 1

' Unit 4: Region 1 Fuel Assembly Hex
unit 4
' Cylinder to contain array of fuel channel elements
  cylinder 1 12.5 415.0 0.0
' Graphite hexprism
  hexprism 2 13.076 415.0 0.0
  array 1 1 place 5 5 1 0.0 0.0 0.0
  media 4 1 2 -1
  boundary 2

' ****************************************************************************
'                        REGION 2 FUEL ASSEMBLIES
' ****************************************************************************

' ----- Medium molten salt fuel channel -----

' Unit 5: Region 2 Fuel Assembly - Single Cylindrical Unit
unit 5
' Single medium fuel channel
  cylinder 1 1.2197 415.0 0.0
' Graphite hexprism to contain fuel channel
  hexprism 2 2.6 415.0 0.0
  media 1 1 1
  media 4 1 -1 2
  boundary 2

' Bare Graphite Hexprism is the same for this region as Region 1 (Unit 3)

' Unit 6 - Region 2 Fuel Assembly Hex
unit 6
' Cylinder to contain array of fuel channel elements
  cylinder 1 12.5 415.0 0.0
' Graphite hexprism
  hexprism 2 13.076 415.0 0.0
  array 2 1 place 5 5 1 0.0 0.0 0.0
  media 4 1 2 -1
  boundary 2

' ****************************************************************************
'                        REGION 3 FUEL ASSEMBLIES
' ****************************************************************************

' ----- Large molten salt fuel channel -----

' Unit 7: Region 3 Fuel Assembly - Single Cylindrical Unit
unit 7
' Single large fuel channel
  cylinder 1 1.3485 415.0 0.0
' Graphite hexprism to contain fuel channel
  hexprism 2 2.6 415.0 0.0
  media 1 1 1
  media 4 1 -1 2
  boundary 2

' Bare Graphite Hexprism is the same for this region as Region 1 (Unit 3)

' Unit 8 - Region 3 Fuel Assembly Hex
unit 8
' Cylinder to contain array of fuel channel elements
  cylinder 1 12.5 415.0 0.0
' Graphite hexprism
  hexprism 2 13.076 415.0 0.0
  array 3 1 place 5 5 1 0.0 0.0 0.0
  media 4 1 2 -1
  boundary 2

' ****************************************************************************
'                        GRAPHITE REFLECTOR CELLS
' ****************************************************************************

' Single hexprism of graphite
unit 9
  hexprism 1 13.076 415.0 0.0
  media 4 1 1
  boundary 1

' ****************************************************************************
'                              GLOBAL UNIT
' ****************************************************************************

global unit 10
' Core Volume (fill with Array 4)
  cylinder 1 190 415 0.0
' Hastelloy N Core Blanket
  cylinder 2 195 415 0.0
' Downcomer Region
  cylinder 3 200 440.5 -10.5
' Hastelloy N Reactor Vessel
  cylinder 4 205 565 -15.5
' Reflector (top) (Hastelloy N)
  cylinder 5 195 430 415
' Outlet Plenum
  cylinder 6 195 435.5 430
' Inlet Plenum
  cylinder 7 195 0.0 -5.5
' Fuel Inside Helium Cylinder
  cylinder 8 200 {h} 440.5
' Helium Cylinder
  cylinder 9 200 560 {h}
' Main array placement
  array 4 1 place 10 10 1 0.0 0.0 0.0
  media 2 1 -1 2
' Downcomer media
  media 1 1 3 -2 -5 -6 -7 -8 -9
' Reactor Vessel Media
  media 2 1 4 -5 -3 -6 -7 -8 -9
' Reflector (top) media
  media 2 1 5
' Outlet Plenum media
  media 1 1 6
' Inlet Plenum media
  media 1 1 7
' Fuel Inside Helium Cylinder Media
  media 1 1 8
' Helium Cylinder media
  media 6 1 9
  boundary 4
end geometry

read array
' Array 1: Fuel Region 1 Single Hex
  ara=1
  prt=yes
  typ=triangular
  nux=9
  nuy=9
  nuz=1
  fill
' z = 1
    3 3 3 3 3 3 3 3 3
     3 3 3 3 3 3 3 3 3
      3 3 3 3 2 2 2 3 3
       3 3 3 2 2 2 2 3 3
        3 3 2 2 2 2 2 3 3
         3 3 2 2 2 2 3 3 3
          3 3 2 2 2 3 3 3 3
           3 3 3 3 3 3 3 3 3
            3 3 3 3 3 3 3 3 3
  end fill

' Array 2: Fuel Region 2 Single Hex
  ara=2
  prt=yes
  typ=triangular
  nux=9
  nuy=9
  nuz=1
  fill
' z = 1
    3 3 3 3 3 3 3 3 3
     3 3 3 3 3 3 3 3 3
      3 3 3 3 5 5 5 3 3
       3 3 3 5 5 5 5 3 3
        3 3 5 5 5 5 5 3 3
         3 3 5 5 5 5 3 3 3
          3 3 5 5 5 3 3 3 3
           3 3 3 3 3 3 3 3 3
            3 3 3 3 3 3 3 3 3
  end fill

' Array 3: Fuel Region 2 Single Hex
  ara=3
  prt=yes
  typ=triangular
  nux=9
  nuy=9
  nuz=1
  fill
' z = 1
    3 3 3 3 3 3 3 3 3
     3 3 3 3 3 3 3 3 3
      3 3 3 3 7 7 7 3 3
       3 3 3 7 7 7 7 3 3
        3 3 7 7 7 7 7 3 3
         3 3 7 7 7 7 3 3 3
          3 3 7 7 7 3 3 3 3
           3 3 3 3 3 3 3 3 3
            3 3 3 3 3 3 3 3 3
  end fill

' Array 4: Main Array Containing All Hex Elements
  ara=4
  prt=yes
  typ=triangular
  nux=21
  nuy=21
  nuz=1
  fill
' z = 1
    9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
     9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
      9 9 9 9 9 9 9 9 9 9 9 9 8 8 9 9 9 9 9 9 9
       9 9 9 9 9 9 9 9 9 8 8 8 8 8 8 8 9 9 9 9 9
        9 9 9 9 9 9 9 9 8 6 6 6 6 6 6 8 9 9 9 9 9
         9 9 9 9 9 9 8 8 6 6 6 6 6 6 6 8 9 9 9 9 9
          9 9 9 9 9 8 8 6 6 6 6 6 6 6 6 8 8 9 9 9 9
           9 9 9 9 9 8 6 6 6 4 1 4 6 6 6 8 8 9 9 9 9
            9 9 9 9 8 6 6 6 1 4 4 1 6 6 6 8 9 9 9 9 9
             9 9 9 8 6 6 6 4 4 4 4 4 6 6 6 8 9 9 9 9 9
              9 9 9 8 6 6 6 1 4 4 1 6 6 6 8 9 9 9 9 9 9
               9 9 9 8 6 6 6 4 1 4 6 6 6 8 9 9 9 9 9 9 9
                9 9 8 8 6 6 6 6 6 6 6 6 8 8 9 9 9 9 9 9 9
                 9 9 8 8 6 6 6 6 6 6 6 8 8 9 9 9 9 9 9 9 9
                  9 9 9 8 6 6 6 6 6 6 8 9 9 9 9 9 9 9 9 9 9
                   9 9 9 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9
                    9 9 9 9 9 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9
                     9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
                      9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
                       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
                        9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
  end fill
end array

end data
end
        '''

            fout = open(self.deck_name, "w")  # Dump deck into file
            fout.write(keno_deck)
            fout.close()
            os.chdir('..')


    def write_qsub_file(self):
        'Writes a single qsub file outside of the refuel directories.'
        qsub_content = '''
#!/bin/bash

#PBS -V
#PBS -q fill
#PBS -l nodes=1:ppn=64

cd $PBS_O_WORKDIR

hostname

module unload mpi
module load openmpi/2.1.6
#module load scale/6.3.b15
#module load scale/ondrejch
module load scale
export DATA=/opt/scale6.3_data

#export scale_input_dir=/home/ondrejch/projects/myMSR/SCALE
export HDF5_USE_FILE_LOCKING=FALSE

scalerte -m -N $PBS_NUM_PPN EIRENE.inp
grep 'best estimate' refuel_*/EIRENE.out | sed -E -e s/^refuel.//g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g  > data-temp.out'''
        try:                # Write the deck
            f = open(self.qsub_path, 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to file", f)
            print(e)

    def run_deck(self):
        'Runs the KENO decks in parallel.'
        for x in range(len(self.rvols)):
            deckpath = f'refuel_{self.rvols[x]:5.01f}'
            os.chdir(deckpath)
            os.system('qsub ~/EIRENE/Scripts/Sourdough/runEIRENE-Scale.sh')
            os.chdir('..')

    def read_outfile(self):
        'Checks whether the run is finished and reads k-eff data.'
        if os.path.exists(self.deck_path+'/data-temp.out'):
          pass
        else:                 
            return False
        filename = "data-temp.out"
        data = genfromtxt("{}".format(filename), delimiter='')
        return data
    
    def get_crit_refuel(self):
        'Fits the k-eff data and returns the critical refuel amount.'
        filename = "data-temp.out"
        data = self.read_outfile(filename)
        refuel = data[:,0]     # Refuel amounts in cm3
        keff = data[:,1]       # k-eff data
        kerr = data[:,2]       # k-eff error
        # Fit the data:
        m, b = np.polyfit(refuel, keff, 1)
        self.crit_refuel = (1 - b)/m     # Refuel amount for which k-eff = 1
        return self.crit_refuel
    
    def plot_refuel_keff(self):
        'Optional: plots the k-eff data for each refuel amount.'
        filename = "data-temp.out"
        data = self.read_outfile(filename)
        refuel = data[:,0]     # Refuel amounts in cm3
        keff = data[:,1]       # k-eff data
        kerr = data[:,2]       # k-eff error
        # Fit the data:
        m, b = np.polyfit(refuel, keff, 1)
        fit = m*refuel + b
        plt.plot(refuel, keff, ls='', c='b', marker='o', markersize=4)
        plt.errorbar(refuel, keff, kerr, ls='', c='b', capsize=1)
        plt.plot(refuel, m*refuel + b, c='b')
        plt.grid()
        plt.xlabel("Refuel amount [$cm^{3}$]", fontsize=12)
        plt.ylabel("$k_{eff}$", fontsize=12)
        plt.tight_layout()
        plt.savefig("refuel-keff.png")

#######################################################################################################################################

if __name__ == '__main__':
    print("This module generates and runs parallel KENO decks for different refuel volume additions.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    decks = Refuel_KENO_deck()
    print("***** Writing SCALE KENO decks...")
    decks.write_KENO_decks()
    decks.write_qsub_file()
    print("All done! Submitting jobs...")
    decks.run_decks()
    while not decks.read_outfile():
        print("The atoms are still working; please stand by...")
        time.sleep(20.0)
    decks.get_crit_refuel()
    decks.plot_refuel_keff()
    print(decks.crit_refuel)

##########################################################################################################################################