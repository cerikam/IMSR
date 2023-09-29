# ***************************************************************************************************
#
#                            EIRENE Sourdough Initialize and Run BOC Core
#
#   By: C. Erika Moss and Dr. Ondrej Chvala
#
# ****************************************************************************************************

import time
import os
import shutil
import salts_wf

# * Initialize BOC TRITON core and submit job to run TRITON.
# * Sleep until TRITON run is finished.
# * Read in .f71 file data from TRITON with write_burned_salt.py. Write new material input deck with write_refuel.py and mix_salts.py at different refuel amounts.
# * Submit KENO jobs in parallel.
# * Sleep until KENO jobs are finished.
# * Submit convert-data.sh shell script to extract best estimate k-eff data from finished KENO jobs.
# * Fit k-eff data and return critical refuel amount with get_crit_refuel.py
# * Write new TRITON deck using critical refuel amount.
# * Repeat process for new depletion step (with the TRITON deck from the previous step being the new BOC core).

#####################################################
#           Initialize BOC TRITON Core
#####################################################

# Initialize and create the input deck

class BOC_core(object):
    'Initial EIRENE TRITON deck.'
    def __init__(self):
        ##################################
        # BOC core and running parameters
        ##################################

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
        #    Set TRITON Params.
        # ************************
        self.dep_step:float = 7          # Length of single depletion step in days
        self.npg:float = 40000           # Number per generation (npg)
        self.nsk:float = 50              # Number of skipped generations (nsk)
        self.gen:float = 5000            # Number of generations (gen)
        self.sig:float = 50e-5           # (sig)
        
        # ************************
        #    Running Parameters
        # ************************

        self.queue:str = 'fill'                     # NECluster queue
        self.ppn:int = 64                           # ppn core count
        self.deck_name:str = 'EIRENE.inp'           # TRITON input file name
        self.deck_path:str = '.'                    # Where to run the TRITON deck
        self.qsub_path:str = os.path.expanduser('~/EIRENE/Scripts/Sourdough/runEIRENE-Scale.sh')  # Full path to the qsub script

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
        
    def write_scale_mat(self):
        s = salts_wf.Salt(self.saltmix(self.UF4molpct), self.Uenrichfrac)
        scale_fuel = s.scale_mat(self.fs_tempK, self.fs_tempK, 1, 1)
        return scale_fuel

    def write_deck(self) -> str:
        'Writes a SCALE input deck for the initial BOC core'
        scale_fuel = self.write_scale_mat()
        triton_deck = f'''=t6-depl parm=(addnux=4)
EIRENE SCALE/TRITON model, UF4 mol% = {self.UF4molpct}
ce_v7.1

read comp
{scale_fuel}

' Hastelloy N
wtptHastelloy 2 {self.hast_dens} 5
        28000 71.0
        24000 7.0
        26000 5.0
        14000 1.0
        42000 16.0
        1.0 {self.mat_tempK}
        28058 67.6 28060 32.4
        24052 100.0
        26056 100.0
        14028 100.0
        42092 14.65 42094 9.19 42095 15.87 42096 16.67 42097 9.58 42098 24.29 42100 9.75 end
' Gd Shutdown Rods
gd 3 den=7.901000 1.0 923.15
    64158 14.8
    64156 10.5
    64155 4.8
    64160 11.9
    64157 55.7
    64152 2.3 end
' Graphite
graphite 4 den={self.gr_dens} 1.0 {self.gr_tempK} end
' Stainless Steel SS316
ss316 5 den={self.ss316_dens} 1.0 {self.mat_tempK} end
' Helium gas
he 6 den={self.he_dens} 1.0 {self.mat_tempK} end
'Dummy composition for off gas removal
xe-135     11    0    1.00000E-20    300    end
kr-85      11    0    1.00000E-20    300    end
'Dummy composition for noble metals
xe-135     12    0    1.00000E-20    300    end
kr-85      12    0    1.00000E-20    300    end
end comp

' Mixtures to deplete
read depletion
1 decayonly 11 12 end
end depletion

' Burnup values. Power in MW/MTiHM = 400.0 / 6.883748384818199 = 58.10787635442664
read burndata
power=58.10787635442664 burn={self.dep_step} nlib=1 end
end burndata

read timetable
'removal of gaseous fission products
'removal time: 50 s
'removal rate: 1/(50s)
flow
from 1 to 11
type fractional_removal
units pers
nuclides   Kr          Xe          end
rates      2.00000E-02 2.00000E-02 end
time       0.0 end
multiplier 1.0 end
end flow
'removal of gaseous fission products
'removal time: 2.4 h
'removal rate: 1/(2.4h * 3600s)
flow
from 1 to 12
type fractional_removal
units pers
nuclides   Zn          Ga          Ge          As
            Nb          Mo          Tc          Ru
            Rh          Pd          Ag          Cd
            In          Sn          Sb          end
rates      1.15741E-04 1.15741E-04 1.15741E-04 1.15741E-04
            1.15741E-04 1.15741E-04 1.15741E-04 1.15741E-04
            1.15741E-04 1.15741E-04 1.15741E-04 1.15741E-04
            1.15741E-04 1.15741E-04 1.15741E-04 end
time       0.0 end
multiplier 1.0 end
end flow
end timetable

read model

read parameters
npg={self.npg} nsk={self.nsk} gen={self.gen} sig={self.sig}
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
' Helium Cylinder
cylinder 8 200 560 440.5
' Main array placement
array 4 1 place 10 10 1 0.0 0.0 0.0
media 2 1 -1 2
' Downcomer media
media 1 1 3 -2 -5 -6 -7 -8
' Reactor Vessel Media
media 2 1 4 -5 -3 -6 -7 -8
' Reflector (top) media
media 2 1 5
' Outlet Plenum media
media 1 1 6
' Inlet Plenum media
media 1 1 7
' Helium Cylinder media
media 6 1 8
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

end model

end
'''
        return triton_deck
    
  
    def save_deck(self):
        'Saves the SCALE TRITON deck as a .inp file.'
        try:
            os.makedirs(self.deck_path, exist_ok=True)
            fh = open(self.deck_path + '/' + self.deck_name, 'w')
            fh.write(self.write_deck())
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ",
                  self.deck_path + '/' + self.deck_name)
            print(e)
    
    def write_qsub_file(self):
        'Writes a submission script for submitting SCALE TRITON job.'
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

scalerte -m -N $PBS_NUM_PPN EIRENE.inp'''
        try:                # Write the deck
            f = open(self.qsub_path, 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to file", f)
            print(e)
        
    def run_deck(self):
        'Runs the deck using qsub_path script'
        if self.queue == 'local':    # Run the deck locally
            os.chdir(self.deck_path)
            os.system(self.qsub_path)
        else:               # Submit the job on the cluster
            os.system('cd ' + self.deck_path + ' && qsub ' + self.qsub_path)

if __name__ == '__main__':
    print("This module generates and runs a BOC EIRENE core for TRITON.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    boc = BOC_core()
    print("***** SCALE TRITON deck: \n" + boc.write_deck() + "\n***** ")
#    mylat.deck_path = os.path.expanduser('~/tmp/lat_test')
    boc.write_qsub_file()
    boc.save_deck()
    boc.run_deck()
#    while not mylat.get_calculated_values():
#        print("Wating for Serpent ...")
#        time.sleep(5.0)
#    print(mylat.k, mylat.cr)