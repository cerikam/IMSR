# ******************************************************************************************************
#
#                               EIRENE Get Critical Refuel Amount
#
#   Purpose: Read k-eff values from KENO runs and fit for criticality. Return refuel amount which
#            results in k-eff = 1
#
# *******************************************************************************************************

# FIRST: Submit jobs with qsub line in write_KENO_deck.py
# SECOND: Wait for jobs to finish. Once finished, have another qsub to submit shell script with "grep best estimate" and creates
#         a data-temp.out file
# THIRD: Read results of data-temp.out and plot data
# FOURTH: Fit data, get refuel amount that makes k-eff = 1
# FIFTH: Write new TRITON deck using refuel amount from prev.

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

# First: Have function that does qsub convert-data.sh

def read_outfile(filename):
    '''Reads the data-temp.out file created from running the shell script.'''
    data = genfromtxt("{}".format(filename), delimiter='')
    return data

# Sort, plot, and fit k-eff data:

def fit_refuel_keff(filename):
    '''Sorts, plots, and fits k-eff data with respect to refuel amount in cm3.
        Data is read using read_outfile function.'''
    data = read_outfile(filename)
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
    crit_refuel = (1 - b)/m     # Refuel amount for which k-eff = 1
    return crit_refuel

test = fit_refuel_keff("data-temp.out")
print(test)
