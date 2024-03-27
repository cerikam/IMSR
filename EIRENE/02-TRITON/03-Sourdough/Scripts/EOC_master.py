###########################################################################################
#
#          EOC Master Script: Create and Run EOC Decks for Each Depletion Step
#
#     By: C. Erika Moss
#
#     Purpose: Writes and runs EOC KENO jobs for EIRENE (this script is intended to be
#              used with EOC_deck.py on the NE Cluster rather than on a local machine)
#
###########################################################################################

import time
import os
import shutil
import salts_wf
import initialize_BOC
from initialize_BOC import BOC_core
import numpy as np
from numpy import genfromtxt
import math
import sys, re
import subprocess
import EOC_deck
from EOC_deck import EOC_Deck

if __name__ == '__main__':
    print("This master EOC script creates and runs an EOC deck for each depletion step (BOC excluded).")
    input("Press Ctrl+C to quit, or press enter to run it.")
    eoc = EOC_Deck()
    iter = np.arange(1, 301, 1)        # Depletion steps
    for i in iter:
        run_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(i))
        os.chdir(run_path)
        eoc.write_EOC_KENO(i)
        eoc.write_qsub_file(i)
        eoc.run_SCALE(i)
    print("*******************************************")
    print("All done! The atoms are very happy now :D")
    print("*******************************************")
