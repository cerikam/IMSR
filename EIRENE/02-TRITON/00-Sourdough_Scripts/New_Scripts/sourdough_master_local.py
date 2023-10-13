# ****************************************************************************************************
#
#                             EIRENE Sourdough Master Script (Local)
#
#    By C. Erika Moss and Dr. Ondrej Chvala
#
# *****************************************************************************************************

import initialize_BOC_local
from initialize_BOC_local import BOC_core
import run_KENO_locally
from run_KENO_locally import Refuel_Deck
import time
import os
import shutil
import salts_wf
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
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

    boc = BOC_core()
    print("Writing BOC EIRENE TRITON deck...")
    boc.save_deck()
    boc.write_sub_shell()
    boc.add_shell_permission()
    print("Running TRITON deck...")
    boc.run_SCALE()

    ###################################################################
    #                  Cycle Through Depletion Steps 
    ###################################################################

    decks = Refuel_Deck()
    iter = np.linspace(1,3,3)        # Depletion steps
    for i in iter:
        decks.write_KENO_decks(i)
        decks.write_sub_shell(i)
        decks.write_conv_data(i)
        decks.add_shell_permission(i)
        decks.add_shell_permission_cd(i)
        print("All done! Running job for depletion step {}".format(i))
        decks.run_SCALE_KENO(i)
        print("Finished running KENO for depletion step {}.".format(i))
        print("Extracting k-eff data for depletion step {}...".format(i))
        decks.convert_data(i)
        print("Data written to data-temp.out file.")
        decks.read_outfile(i)
        decks.get_crit_refuel(i)
        decks.plot_refuel_keff(i)
        print("Critical refuel amount for depletion step {}:".format(i))
        print(decks.crit_refuel)
        print("Writing new TRITON deck with critical refuel amount...")
        decks.write_new_TRITON_deck(i)
        decks.save_deck(i)
        print("***** SCALE TRITON deck: \n" + decks.write_new_TRITON_deck(i) + "\n***** ")
        print("Running TRITON deck...")
        decks.run_SCALE(i)
    print("*******************************************")
    print("All done! The atoms are very happy now :D")
    print("*******************************************")

#################################################################################################