#!/usr/bin/env python3
#
# Serpent2, \alpha_T
#
# Ondrej Chvala, ochvala@utexas.edu
# MIT license
#

import serpentTools
import os
import numpy as np
import re
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import json

def fitf(x, a, b) -> float:
    return a*x + b


fout = open("kdata.dat","w")

enrs = []
ks = [] ; kerrs = []
rhos = [] ; rhoerrs = []
k0:float = -1.0
k0err:float = -1.0
for enrpct in np.linspace(2.7, 2.9, 11):
    deckpath = f'enr_{enrpct:5.03f}'
    res = serpentTools.read(f'{deckpath}/SIRENE_res.m')
    k:float    = float(res.resdata['absKeff'][0])
    kerr:float = float(res.resdata['absKeff'][1])*k
    rho:float    = 1e5 * (k - 1.0) / k
    rhoerr:float = float(res.resdata['absKeff'][1]) * rho
    fout.write(f'{enrpct:5.03f}   {k:7.5f} {kerr:7.5f}    {rho:.2f} {rhoerr:.2f}\n')
    if (False):         # Six factor formula
        k0 = k
        k0err = kerr
        sixFfEpsilon:float    = float(res.resdata['sixFfEpsilon'][0])
        sixFfEpsilonerr:float = float(res.resdata['sixFfEpsilon'][1])
        sixFfF:float    = float(res.resdata['sixFfF'][0])
        sixFfFerr:float = float(res.resdata['sixFfF'][1])
        sixFfP:float    = float(res.resdata['sixFfP'][0])
        sixFfPerr:float = float(res.resdata['sixFfP'][1])
        sixFfEta:float    = float(res.resdata['sixFfEta'][0])
        sixFfEtaerr:float = float(res.resdata['sixFfEta'][1])
    enrs.append(enrpct)
    ks.append(k) ; kerrs.append(kerr)
    rhos.append(rho) ; rhoerrs.append(rhoerr)


# Plot
#if intrusion > 0:
#    continue
plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.figsize"] = [7, 7]
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16

# k_eff

# Fit
popt, pcov = curve_fit(fitf, enrs, ks, sigma=kerrs, absolute_sigma=True, method='lm')
perr = np.sqrt(np.diag(pcov))

fig, ax = plt.subplots(1,1)
ax.set_title(f"EIRENE FLiBe-U 5mol% UF$_4$", fontsize=20)
ax.set_xlabel("Uranium enrichment %", fontsize=16)
ax.set_ylabel("SIRENE k$_{eff}$", fontsize=16)
data = ax.errorbar(enrs, ks, xerr=None, yerr=kerrs, fmt='o', color='blue', markersize=3.6)
xfit = np.linspace(enrs[0],enrs[-1],100)
dfit = ax.plot(xfit, fitf(xfit, popt[0], popt[1]), color="green", linewidth=0.4,
       label=f"{popt[0]:.4f} +- {perr[0]:.4f} dK / Uenr%")
leg = ax.legend()
plt.grid()
#plt.show()
fig.savefig(f"plot_enr_k.png", bbox_inches="tight", facecolor='white')
plt.close()

# rho [pcm]

# Fit
popt, pcov = curve_fit(fitf, enrs, rhos, sigma=rhoerrs, absolute_sigma=True, method='lm')
perr = np.sqrt(np.diag(pcov))

fig, ax = plt.subplots(1,1)
ax.set_title(f"EIRENE FLiBe-U 5mol% UF$_4$", fontsize=20)
ax.set_xlabel("Uranium enrichment %", fontsize=16)
ax.set_ylabel("SIRENE ρ [pcm]", fontsize=16)
data = ax.errorbar(enrs, rhos, xerr=None, yerr=rhoerrs, fmt='o', color='red', markersize=3.6)
xfit = np.linspace(enrs[0],enrs[-1],100)
dfit = ax.plot(xfit, fitf(xfit, popt[0], popt[1]), color="green", linewidth=0.4,
        label=f"{popt[0]:.2f} +- {perr[0]:.2f} pcm / Uenr%")
leg = ax.legend()
plt.grid()
#plt.show()
fig.savefig(f"plot_enr-rho.png", bbox_inches="tight", facecolor='white')
plt.close()


