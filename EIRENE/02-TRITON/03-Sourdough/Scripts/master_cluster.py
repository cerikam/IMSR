# ****************************************************************************************************
#
#                             EIRENE Sourdough Master Script (Local)
#
#    By C. Erika Moss and Dr. Ondrej Chvala
#
# *****************************************************************************************************

import initialize_BOC_local
from initialize_BOC_local import BOC_core
import run_sourdough
from run_sourdough import Refuel_Deck
import time
import os
import shutil
import salts_wf
import numpy as np
from numpy import genfromtxt
#import matplotlib.pyplot as plt
import math
import sys, re
import subprocess

#################################################################################################

if __name__ == '__main__':
    print("This master script runs a full sourdough depletion simulation for EIRENE.")
    input("Press Ctrl+C to quit, or enter else to test it.")

    ##################################################################
    #        Initialize the BOC Core and Run the TRITON Deck
    ##################################################################
    
    # NOTE: Uncomment the lines in this section to also run a BOC deck.

#    boc = BOC_core()
#    print("Writing BOC EIRENE TRITON deck...")
#    boc.save_deck()
#    boc.write_sub_shell()
#    boc.add_shell_permission()
#    print("Running TRITON deck...")
#    boc.run_SCALE()

    ###################################################################
    #                  Cycle Through Depletion Steps 
    ###################################################################

    decks = Refuel_Deck()
    iter = np.arange(1, 301, 1)        # Depletion steps
    for i in iter:
        decks.write_KENO_decks(i)
        decks.write_conv_data(i)
        decks.add_shell_permission_cd(i)
        time.sleep(1700)
        print("Extracting k-eff data for depletion step {}...".format(i))
        decks.convert_data(i)
        print("Data written to data-temp.out file.")
        decks.read_outfile(i)
        decks.get_crit_refuel(i)
        print("Critical refuel amount for depletion step {}:".format(i))
        print(decks.crit_refuel)
        print("Writing new TRITON deck with critical refuel amount...")
        decks.write_new_TRITON_deck(i)
        decks.write_qsub_file(i)
        print("Running TRITON deck...")
        decks.run_SCALE(i)
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(i))
        os.chdir(main_path)
        check = False       # Initialize check
        while check == False:
          if os.path.exists('./EIRENE.f71'):
            check = True
            break
          else:                 
              check = False
              print("The atoms are still working; please stand by...")
              time.sleep(10.0)
    print("*******************************************")
    print("All done! The atoms are very happy now :D")
    print("*******************************************")

#################################################################################################
