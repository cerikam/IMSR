#!/bin/env python3


""" Serpent model of example Integral REactor for Nuclear Education (SIRENE)
Script to find Fuel Temerature Coefficient (FTC) of Reactivity
Ondrej Chvala <ochvala@utk.edu>
MIT license
"""

import salts
import os
import numpy as np
import h5py

def s2comp(h5file:str="") -> str:
    ''' Reads HDF5 SHIFT file and returns material deck for Serpent.
    Only the relevant ones are parsed: fuel salt and moderator graphite.
    '''
    out:str = ''
    Shift2SerpentComps = {
    1 : ['fuel'],
    2 : ['hastelloy'],
    3 : ['gd'],
    4 : ['graphite'],
#    5 : ['steel'],
    6 : ['helium'],
    }
    f = h5py.File(h5file)
    comp = f['comp']['compositions']

    for i in range(1, comp.len()):
        if not i in Shift2SerpentComps.keys():
            continue
        c = comp[i]
        T = c[2]
        dens = c[3]
        isos = c[4]
        s2isomat:str = ''
        #print(f"*** {i}, temp {T} K, rho {dens} cc")
        grmassfrac:float = 0;
        if T > 900:
            xlib_fuel = '82c'
        if T < 900:
            xlib_fuel = '81c'
        if T < 600:
            xlib_fuel = '80c'
        for iso in isos:            # Loop over isotopes in mixture and make Serpent-like entries
            if iso[1] <= 0.0:
                continue
            if int(iso[0]) > 3006000:   # Special treatment for graphite
                grmassfrac += float(iso[1])
                continue
            isomat = f"{iso[0]:-7d}.{xlib_fuel} -{iso[1]:021.19f}"
            #print(isomat)
            s2isomat = s2isomat + f"{isomat}\n"
        if grmassfrac > 0.0:        # There was graphite!
            isomat = f"{6000:-7d}.{xlib} -{grmassfrac:021.19f}"
            #print(isomat)
            s2isomat = s2isomat + f"{isomat}\n"
        for m in Shift2SerpentComps[i]:  # Loop ever Serpent materials to generate
            #print(m)
            s2mathead:str = ''
            if grmassfrac > 0.0:        # There was graphite!
                s2mathead = f"mat {m}  -{dens} tmp {T} moder grph 6000"
            else:
                s2mathead = f"mat {m}  -{dens} tmp {T}"
            out = out + f"{s2mathead}\n{s2isomat}\n"
    return out


T   = 923.15          # material temperature [K]
Tgr = T


for salt_tempC in np.linspace(600, 700, 11):
    deckpath = f'FTC_{salt_tempC:5.01f}'
    if not os.path.isdir(deckpath):
        os.mkdir(deckpath)
    os.chdir(deckpath)

    xlib = '82c'    # material temperature [K]
    if T < 900:
        xlib = '81c'
    if T < 600:
        xlib = '80c'
    

    Tgr = T # Sab library selection based on graphite T
    grlib0 = '16t'
    grlib1 = '17t'
    if Tgr<=1000.0:
        grlib0 = '15t'
        grlib1 = '16t'
    if Tgr<=800.0:
        grlib0 = '14t'
        grlib1 = '15t'
    if Tgr<=700.0:
        grlib0 = '13t'
        grlib1 = '14t'
    if Tgr<=600.0:
        grlib0 = '12t'
        grlib1 = '13t'
    if Tgr<=500.0:
        grlib0 = '11t'
        grlib1 = '12t'
    if Tgr<=400.0:
        grlib0 = '10t'
        grlib1 = '11t'



    s2matdeck:str = s2comp(f'../../../EIRENE/10-SHIFT/02-FTC/{deckpath}/EIRENE.shift-output.h5')

    fout = open('SIRENE','w')
    fout.write(f'''
% *********************************************************************
%
%
%                EIRENE REACTOR SERPENT-2 MODEL
%
%
% *********************************************************************


% *********************************************************************
%                    REGION 1 FUEL ASSEMBLIES
% *********************************************************************

% -------------------------------------------------------
%              Small molten salt fuel channels
% -------------------------------------------------------

% --- Universe 1: Single small fuel channel in graphite hex ---

% Surfaces:

surf 1 cyl 0 0 1.157 0 415       % Single small fuel channel
surf 2 hexxprism 0 0 2.6 0 415   % Graphite hex to contain channel

% Cells:

cell 1 1 fuel -1
cell 2 1 graphite 1 -2
cell 3 1 graphite 2

% --- Universe 2: Bare graphite hexprism for creating the array

surf 3 hexxprism 0 0 2.6 0 415   % Bare grahpite hex

cell 4 2 graphite -3
cell 5 2 graphite 3

lat R1 2 0.0 0.0 9 9 5.2
2 2 2 2 2 2 2 2 2
 2 2 2 2 2 2 2 2 2
  2 2 2 2 1 1 1 2 2
   2 2 2 1 1 1 1 2 2
    2 2 1 1 1 1 1 2 2
     2 2 1 1 1 1 2 2 2
      2 2 1 1 1 2 2 2 2
       2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2

% --- Universe S: Single Small Fuel Assembly Hex ---

surf 4 hexxprism 0 0 13.076 0 415  % Hex to contain R1 channels

cell 6 S fill R1 -4
cell 7 S graphite 4

% *********************************************************************
%                    REGION 2 FUEL ASSEMBLIES
% *********************************************************************

% -------------------------------------------------------
%              Medium molten salt fuel channels
% -------------------------------------------------------

% --- Universe 3: Single medium fuel channel in graphite hex ---

surf 5 cyl 0 0 1.2197 0 415       % Single medium fuel channel
surf 6 hexxprism 0 0 2.6 0 415    % Graphite hex to contain channel

cell 8 3 fuel -5
cell 9 3 graphite 5 -6
cell 10 3 graphite 6

lat R2 2 0.0 0.0 9 9 5.2
2 2 2 2 2 2 2 2 2
 2 2 2 2 2 2 2 2 2
  2 2 2 2 3 3 3 2 2
   2 2 2 3 3 3 3 2 2
    2 2 3 3 3 3 3 2 2
     2 2 3 3 3 3 2 2 2
      2 2 3 3 3 2 2 2 2
       2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2

% --- Universe M: Single Medium Fuel Assembly Hex

surf 7 hexxprism 0 0 13.076 0 415  % Hex to contain R2 channels

cell 11 M fill R2 -7
cell 12 M graphite 7

% *********************************************************************
%                    REGION 3 FUEL ASSEMBLIES
% *********************************************************************

% -------------------------------------------------------
%              Large molten salt fuel channels
% -------------------------------------------------------

% --- Universe 4: Single large fuel channel in graphite hex ---

surf 8 cyl 0 0 1.3485 0 415       % Single large fuel channel
surf 9 hexxprism 0 0 2.6 0 415    % Graphite hex to contain channel

cell 13 4 fuel -8
cell 14 4 graphite 8 -9
cell 15 4 graphite 9

lat R3 2 0.0 0.0 9 9 5.2
2 2 2 2 2 2 2 2 2
 2 2 2 2 2 2 2 2 2
  2 2 2 2 4 4 4 2 2
   2 2 2 4 4 4 4 2 2
    2 2 4 4 4 4 4 2 2
     2 2 4 4 4 4 2 2 2
      2 2 4 4 4 2 2 2 2
       2 2 2 2 2 2 2 2 2
        2 2 2 2 2 2 2 2 2

% --- Universe L: Single Large Fuel Assembly Hex

surf 10 hexxprism 0 0 13.076 0 415  % Hex to contain R3 channels

cell 16 L fill R3 -10
cell 17 L graphite 10

% *********************************************************************
%                      BARE GRAPHITE HEX CELL
% *********************************************************************

% Universe G - Bare Graphite Hex

surf 11 hexxprism 0 0 13.076 0 415

cell 18 G graphite -11
cell 19 G graphite 11

% *********************************************************************
%                          CONTROL RODS
% *********************************************************************

% Universe C - Control Rods

surf 100 cyl 0 0 7.0 130 415           % Shutdown rod (Gd)
surf 101 cyl 0 0 7.2 130 415           % Guide tube
surf 102 cyl 0 0 7.2 0 130             % Fuel region below
surf 103 hexxprism 0 0 13.076 0 415    % Graphite hexprism to contain rod

cell 100 C gd -100                     % Shutdown rod (Gd)
cell 101 C hastelloy 100 -101          % Guide tube (Hastelloy)
cell 102 C fuel -102                   % Fuel region below
cell 103 C graphite 100 101 102 -103   % Graphite hexprism to contain rod
cell 104 C graphite 103


% *********************************************************************
%                         MAIN ASSEMBLY
% *********************************************************************

% ------------------------ LATTICE CODE -------------------------------
%
%          G = bare graphite hexprism
%          L = large fuel channel (fuel region 3)
%          M = medium fuel channel (fuel region 2)
%          S = small fuel channel (fuel region 1)
%          C = control rod
%
% ---------------------------------------------------------------------

lat MAIN 2 0.0 0.0 21 21 26.153
G G G G G G G G G G G G G G G G G G G G G
 G G G G G G G G G G G G G G G G G G G G G
  G G G G G G G G G G G G G G G G G G G G G
   G G G G G G G G G G G G G L L G G G G G G
    G G G G G G G G G G L L L L L L L G G G G
     G G G G G G G G G L M M M M M M L G G G G
      G G G G G G G L L M M M M M M M L L G G G
       G G G G G G L L M M M M M M M M L L G G G
        G G G G G G L M M M S C S M M M L G G G G
         G G G G G L M M M C S S C M M M L G G G G
          G G G G L M M M S S S S S M M M L G G G G
           G G G G L M M M C S S C M M M L G G G G G
            G G G G L M M M S C S M M M L G G G G G G
             G G G L L M M M M M M M M L L G G G G G G
              G G G L L M M M M M M M L L G G G G G G G
               G G G G L M M M M M M L G G G G G G G G G
                G G G G L L L L L L L G G G G G G G G G G
                 G G G G G G L L G G G G G G G G G G G G G
                  G G G G G G G G G G G G G G G G G G G G G
                   G G G G G G G G G G G G G G G G G G G G G
                    G G G G G G G G G G G G G G G G G G G G G

surf 12 cyl 0 0 190 0 415         % Clinder to contain main array
surf 13 cyl 0 0 195 0 415         % Hastelloy N core blanket
surf 14 cyl 0 0 200 -10.5 440.5   % Downcomer region
surf 15 cyl 0 0 195 415 430       % Hastelloy N reflector
surf 16 cyl 0 0 195 430 435.5     % Outlet plenum
surf 17 cyl 0 0 195 -5.5 0        % Inlet plenum
surf 18 cyl 0 0 200 440.5 560     % Helium cylinder
surf 19 cyl 0 0 205 -15.5 565     % Hastelloy N reactor vessel

cell 20 0 fill MAIN -12                  % Core region (main array fill)
cell 21 0 hastelloy 12 -13               % Hastelloy N core blanket
cell 22 0 fuel 12 13 15 16 17 18 -14     % Downcomer region
cell 23 0 hastelloy 12 13 -15            % Hastelloy N reflector
cell 24 0 fuel 12 13 -16                 % Outlet plenum
cell 25 0 fuel 12 13 -17                 % Inlet plenum
cell 26 0 helium 12 13 -18               % Helium cylinder
cell 27 0 hastelloy -19 15 14 16 17 18   % Hastelloy N reactor vessel
cell 99 0 outside 19

% *********************************************************************
%                             MATERIALS
% *********************************************************************

{s2matdeck}

% Thermal scattering data
therm grph {Tgr} "grph.{grlib0}" "grph.{grlib1}" % Graphite at {Tgr} K


% --- Neutron population and criticality cycles:
set pop 50000 820 20

% --- Geometry and mesh plots:
plot 3 1000 1000
plot 1 1000 1000 0.0 -300 300 -60 600
plot 2 1000 1000 0.0 -300 300 -60 600
mesh 3 1000 1000
mesh 1 1000 1000 0 -300 300 -300 300 -60 600
mesh 2 1000 1000 0 -300 300 -300 300 -60 600

% --- Data Libraries
set acelib "/opt/MCNP6.2/MCNP_DATA/xsdir_mcnp6.2_msr783k.sss"
''')

    os.system('qsub ../../utilities/runSIRENE.sh')
    os.chdir('..')
