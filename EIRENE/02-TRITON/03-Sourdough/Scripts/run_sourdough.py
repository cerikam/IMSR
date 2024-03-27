# ******************************************************************************************************
#
#                            EIRENE Run Sourdough Fuel Cycle Simulation
#
#  By: C. Erika Moss and Dr. Ondrej Chvala
#
# *******************************************************************************************************

import time
import os
import shutil
import salts_wf
import initialize_BOC
from initialize_BOC import BOC_core
import numpy as np
from numpy import genfromtxt
#import matplotlib.pyplot as plt
import math
import sys, re
import subprocess

SCALE_bin_path: str = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')

#####################################################
#     Write KENO Decks Based on Refuel Parameters
#####################################################

class Refuel_Deck(object):
    'Create parallel KENO decks for EIRENE for different refuel volume additions. After obtaining critical refuel amount, writes and runs a new TRITON deck.'
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
        self.renrich = 0.07                                    # Refuel enrichment fraction
        self.rvols = [1000, 30000, 70000, 110000, 160000]      # Refuel volumes in cm^3
        self.V0 = 1.0700515E+07                                # Total fuel salt volume for BOC core in cm^3

        # ************************
        #    Set KENO Params.
        # ************************
        self.npg:float = 20000           # Number per generation (npg)
        self.nsk:float = 20              # Number of skipped generations (nsk)
        self.gen:float = 5020            # Number of generations (gen)
        self.sig:float = 50e-5           # (sig)
        self.dep_step:float = 21          # Length of single depletion step in days

        # ************************
        #    Running Parameters
        # ************************
        self.queue:str = 'fill'                     # NECluster queue
        self.ppn:int = 64                           # ppn core count
        self.deck_name:str = 'EIRENE.inp'           # KENO input file name
        self.f71_name:str = 'EIRENE.f71'            # Name of .f71 TRITON file to read
        self.qsub_name:str = 'runEIRENE-Scale.sh'   # Name of each qsub file
        self.qsub_path:str = os.path.expanduser('~/EIRENE11/SD7pct/runEIRENE-Scale.sh')     # Full path to the qsub script
        self.conv_data_path:str = os.path.expanduser('~/EIRENE11/SD7pct/convert-data.sh')   # Full path to convert-data.sh script
        self.deck_path:str = os.path.expanduser('~/EIRENE11/SD7pct/')                    # Path of main run directory
        self.BOC_f71_path:str = os.path.expanduser('~/EIRENE11/SD7pct/BOC')
        self.power:float = 400.0      # TRITON power in MW (used for calculating power density)
        self.init_MTiHM:float = 6.883864319048779        # Initial MTiHM for BOC core

        ##################################
        #     
        ##################################

    def get_burned_salt_atoms(self, iter):
        'Returns the number of atoms for each constituent of the burned salt.'
        output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=atom", "-idform='{:Ee}{:AAA}{:m}'", self.f71_name], capture_output=True)
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
                densities[nuclide] = float(data[2])  # The [2] here is what causes the code to return only the densities at position 2 of the f71 file (position 0 is isotope names, 1 is BOC data)

        # **** Extract Everything Listed Below: ****
        
        # *** FLiBe Components First for Ease of Identification for Mixing: ***
        
        be9 = densities['be-9']
        f19 = densities['f-19']
        li6 = densities['li-6']
        li7 = densities['li-7']
        u234 = densities['u-234']
        u235 = densities['u-235']
        u236 = densities['u-236']
        u238 = densities['u-238']

        # *********************************

        h1 = densities['h-1']
        h2 = densities['h-2']
        h3 = densities['h-3']
        he3 = densities['he-3']
        he4 = densities['he-4']
        be7 = densities['be-7']
#        be8 = densities['be-8']
        b10 = densities['b-10']
        b11 = densities['b-11']
#        b12 = densities['b-12']
        n14 = densities['n-14']
        n15 = densities['n-15']
#        n16 = densities['n-16']
        o16 = densities['o-16']
        o17 = densities['o-17']
        #f20 = densities['f-20']
        na23 = densities['na-23']
        mg24 = densities['mg-24']
        mg25 = densities['mg-25']
        mg26 = densities['mg-26']
        al27 = densities['al-27']
        si28 = densities['si-28']
        si29 = densities['si-29']
        si30 = densities['si-30']
        p31 = densities['p-31']
        s32 = densities['s-32']
        s33 = densities['s-33']
#        s35 = densities['s-35']
       # s34 = densities['s-34']
        cl35 = densities['cl-35']
        #cl36 = densities['cl-36']
        cl37 = densities['cl-37']
        #cl38 = densities['cl-38']
        ar36 = densities['ar-36']
        #ar37 = densities['ar-37']
        ar38 = densities['ar-38']
        #ar39 = densities['ar-39']
        ar40 = densities['ar-40']
        #ar41 = densities['ar-41']
        #as72 = densities['as-72']
        #as73 = densities['as-73']
        as74 = densities['as-74']
        as75 = densities['as-75']
        #as76 = densities['as-76']
        #as77 = densities['as-77']
        #as78 = densities['as-78']
        #as79 = densities['as-79']
        #as80 = densities['as-80']
        #as81 = densities['as-81']
        #as82 = densities['as-82']
        #as83 = densities['as-83']
        #as84 = densities['as-84']
        #as85 = densities['as-85']
        #as86 = densities['as-86']
        #as87 = densities['as-87']
        #as88 = densities['as-88']
        #as89 = densities['as-89']
        #as90 = densities['as-90']
        k39 = densities['k-39']
        k40 = densities['k-40']
        k41 = densities['k-41']
        ca40 = densities['ca-40']
        #ca41 = densities['ca-41']
        ca42 = densities['ca-42']
        ca43 = densities['ca-43']
        ca44 = densities['ca-44']
        #ca45 = densities['ca-45']
        ca46 = densities['ca-46']
        #ca47 = densities['ca-47']
        ca48 = densities['ca-48']
        #ca49 = densities['ca-49']
        sc45 = densities['sc-45']
        ti46 = densities['ti-46']
        ti47 = densities['ti-47']
        ti48 = densities['ti-48']
        ti49 = densities['ti-49']
        ti50 = densities['ti-50']
        cr50 = densities['cr-50']
        cr52 = densities['cr-52']
        cr53 = densities['cr-53']
        cr54 = densities['cr-54']
        mn55 = densities['mn-55']
        fe54 = densities['fe-54']
        fe56 = densities['fe-56']
        fe57 = densities['fe-57']
        fe58 = densities['fe-58']
        co58 = densities['co-58']
        co59 = densities['co-59']
        ni58 = densities['ni-58']
        ni59 = densities['ni-59']
        ni60 = densities['ni-60']
        ni61 = densities['ni-61']
        ni62 = densities['ni-62']
        ni64 = densities['ni-64']
        cu63 = densities['cu-63']
        #cu64 = densities['cu-64']
        cu65 = densities['cu-65']
        ga69 = densities['ga-69']
        ga71 = densities['ga-71']
        ge70 = densities['ge-70']
       # ge71 = densities['ge-71']
        ge72 = densities['ge-72']
        ge73 = densities['ge-73']
        ge74 = densities['ge-74']
       # ge75 = densities['ge-75']
        ge76 = densities['ge-76']
        se74 = densities['se-74']
        #se75 = densities['se-75']
        se76 = densities['se-76']
        se77 = densities['se-77']
        se78 = densities['se-78']
        se79 = densities['se-79']
        se80 = densities['se-80']
        #se81 = densities['se-81']
        se82 = densities['se-82']
        br79 = densities['br-79']
       # br80 = densities['br-80']
        br81 = densities['br-81']
        kr78 = densities['kr-78']
        #kr79 = densities['kr-79']
        kr80 = densities['kr-80']
        #kr81 = densities['kr-81']
        kr82 = densities['kr-82']
        kr83 = densities['kr-83']
        kr84 = densities['kr-84']
        kr85 = densities['kr-85']
        kr86 = densities['kr-86']
        rb85 = densities['rb-85']
        rb86 = densities['rb-86']
        rb87 = densities['rb-87']
        sr84 = densities['sr-84']
        #sr85 = densities['sr-85']
        sr86 = densities['sr-86']
        sr87 = densities['sr-87']
        sr88 = densities['sr-88']
        sr89 = densities['sr-89']
        sr90 = densities['sr-90']
        y89 = densities['y-89']
        y90 = densities['y-90']
        y91 = densities['y-91']
        #zr89 = densities['zr-89']
        zr90 = densities['zr-90']
        zr91 = densities['zr-91']
        zr92 = densities['zr-92']
        zr93 = densities['zr-93']
        zr94 = densities['zr-94']
        zr95 = densities['zr-95']
        zr96 = densities['zr-96']
        nb93 = densities['nb-93']
        nb94 = densities['nb-94']
        nb95 = densities['nb-95']
        mo92 = densities['mo-92']
        #mo93 = densities['mo-93']
        mo94 = densities['mo-94']
        mo95 = densities['mo-95']
        mo96 = densities['mo-96']
        mo97 = densities['mo-97']
        mo98 = densities['mo-98']
        mo99 = densities['mo-99']
        mo100 = densities['mo-100']
        tc99 = densities['tc-99']
        ru96 = densities['ru-96']
        #ru97 = densities['ru-97']
        ru98 = densities['ru-98']
        ru99 = densities['ru-99']
        ru100 = densities['ru-100']
        ru101 = densities['ru-101']
        ru102 = densities['ru-102']
        ru103 = densities['ru-103']
        ru104 = densities['ru-104']
        ru105 = densities['ru-105']
        ru106 = densities['ru-106']
        rh103 = densities['rh-103']
        #rh104 = densities['rh-104']
        rh105 = densities['rh-105']
        pd102 = densities['pd-102']
        #pd103 = densities['pd-103']
        pd104 = densities['pd-104']
        pd105 = densities['pd-105']
        pd106 = densities['pd-106']
        pd107 = densities['pd-107']
        pd108 = densities['pd-108']
        #pd109 = densities['pd-109']
        pd110 = densities['pd-110']
        ag107 = densities['ag-107']
        #ag108 = densities['ag-108']
        ag109 = densities['ag-109']
        #ag110 = densities['ag-110']
        ag111 = densities['ag-111']
        #cd105 = densities['cd-105']
        cd106 = densities['cd-106']
        #cd107 = densities['cd-107']
        cd108 = densities['cd-108']
        #cd109 = densities['cd-109']
        cd110 = densities['cd-110']
        cd111 = densities['cd-111']
        cd112 = densities['cd-112']
        cd113 = densities['cd-113']
        cd114 = densities['cd-114']
        #cd115 = densities['cd-115']
        cd116 = densities['cd-116']
       # cd117 = densities['cd-117']
       # cd118 = densities['cd-118']
        #cd119 = densities['cd-119']
        #cd120 = densities['cd-120']
        #cd121 = densities['cd-121']
        #cd122 = densities['cd-122']
        #cd123 = densities['cd-123']
        #cd124 = densities['cd-124']
        #cd125 = densities['cd-125']
        #cd126 = densities['cd-126']
        #cd127 = densities['cd-127']
        #cd128 = densities['cd-128']
        #cd129 = densities['cd-129']
        #cd130 = densities['cd-130']
        #cd131 = densities['cd-131']
        #cd132 = densities['cd-132']
        in113 = densities['in-113']
        #in114 = densities['in-114']
        in115 = densities['in-115']
        sn112 = densities['sn-112']
        sn113 = densities['sn-113']
        sn114 = densities['sn-114']
        sn115 = densities['sn-115']
        sn116 = densities['sn-116']
        sn117 = densities['sn-117']
        sn118 = densities['sn-118']
        sn119 = densities['sn-119']
        sn120 = densities['sn-120']
        #sn121 = densities['sn-121']
        sn122 = densities['sn-122']
        sn123 = densities['sn-123']
        sn124 = densities['sn-124']
        sn125 = densities['sn-125']
        sn126 = densities['sn-126']
        sb121 = densities['sb-121']
        #sb122 = densities['sb-122']
        sb123 = densities['sb-123']
        sb124 = densities['sb-124']
        sb125 = densities['sb-125']
        sb126 = densities['sb-126']
        te120 = densities['te-120']
        #te121 = densities['te-121']
        te122 = densities['te-122']
        te123 = densities['te-123']
        te124 = densities['te-124']
        te125 = densities['te-125']
        te126 = densities['te-126']
        #te127 = densities['te-127']
        te128 = densities['te-128']
        #te129 = densities['te-129']
        te130 = densities['te-130']
        #te131 = densities['te-131']
        te132 = densities['te-132']
        i127 = densities['i-127']
        #i128 = densities['i-128']
        i129 = densities['i-129']
        i130 = densities['i-130']
        i131 = densities['i-131']
        #i132 = densities['i-132']
        #i133 = densities['i-133']
        #i134 = densities['i-134']
        i135 = densities['i-135']
        xe123 = densities['xe-123']
        xe124 = densities['xe-124']
        #xe125 = densities['xe-125']
        xe126 = densities['xe-126']
        #xe127 = densities['xe-127']
        xe128 = densities['xe-128']
        xe129 = densities['xe-129']
        xe130 = densities['xe-130']
        xe131 = densities['xe-131']
        xe132 = densities['xe-132']
        xe133 = densities['xe-133']
        xe134 = densities['xe-134']
        xe135 = densities['xe-135']
        xe136 = densities['xe-136']
        #ba128 = densities['ba-128']
        #ba129 = densities['ba-129']
        ba130 = densities['ba-130']
        #ba131 = densities['ba-131']
        ba132 = densities['ba-132']
        ba133 = densities['ba-133']
        ba134 = densities['ba-134']
        ba135 = densities['ba-135']
        ba136 = densities['ba-136']
        ba137 = densities['ba-137']
        ba138 = densities['ba-138']
        #ba139 = densities['ba-139']
        ba140 = densities['ba-140']
        #ba141 = densities['ba-141']
        #ba142 = densities['ba-142']
        #ba143 = densities['ba-143']
        #ba144 = densities['ba-144']
        #ba145 = densities['ba-145']
        #ba146 = densities['ba-146']
        #ba147 = densities['ba-147']
        #ba148 = densities['ba-148']
        #ba149 = densities['ba-149']
        #ba150 = densities['ba-150']
        #ba151 = densities['ba-151']
        #ba152 = densities['ba-152']
        #ba153 = densities['ba-153']
        la138 = densities['la-138']
        la139 = densities['la-139']
        la140 = densities['la-140']
        #ce134 = densities['ce-134']
        #ce135 = densities['ce-135']
        ce136 = densities['ce-136']
        #ce137 = densities['ce-137']
        ce138 = densities['ce-138']
        ce139 = densities['ce-139']
        ce140 = densities['ce-140']
        ce141 = densities['ce-141']
        ce142 = densities['ce-142']
        ce143 = densities['ce-143']
        ce144 = densities['ce-144']
        #ce145 = densities['ce-145']
        #ce146 = densities['ce-146']
        #ce147 = densities['ce-147']
        #ce148 = densities['ce-148']
        #ce149 = densities['ce-149']
        #ce150 = densities['ce-150']
        #ce151 = densities['ce-151']
        #ce152 = densities['ce-152']
        #ce153 = densities['ce-153']
        #ce154 = densities['ce-154']
        #ce155 = densities['ce-155']
        #ce156 = densities['ce-156']
        #ce157 = densities['ce-157']
        pr141 = densities['pr-141']
        pr142 = densities['pr-142']
        pr143 = densities['pr-143']
        nd142 = densities['nd-142']
        nd143 = densities['nd-143']
        nd144 = densities['nd-144']
        nd145 = densities['nd-145']
        nd146 = densities['nd-146']
        nd147 = densities['nd-147']
        sm149 = densities['sm-149']
        sm150 = densities['sm-150']
        sm151 = densities['sm-151']
        sm152 = densities['sm-152']
        sm153 = densities['sm-153']
        sm154 = densities['sm-154']
        eu151 = densities['eu-151']
        eu152 = densities['eu-152']
        eu153 = densities['eu-153']
        eu154 = densities['eu-154']
        eu155 = densities['eu-155']
        eu156 = densities['eu-156']
        eu157 = densities['eu-157']
        gd152 = densities['gd-152']
        gd153 = densities['gd-153']
        gd154 = densities['gd-154']
        gd155 = densities['gd-155']
        gd156 = densities['gd-156']
        gd157 = densities['gd-157']
        gd158 = densities['gd-158']
        #gd159 = densities['gd-159']
        gd160 = densities['gd-160']
        tb159 = densities['tb-159']
        tb160 = densities['tb-160']
        dy156 = densities['dy-156']
      #  dy157 = densities['dy-157']
        dy158 = densities['dy-158']
        #dy159 = densities['dy-159']
        dy160 = densities['dy-160']
        dy161 = densities['dy-161']
        dy162 = densities['dy-162']
        dy163 = densities['dy-163']
        dy164 = densities['dy-164']
        ho165 = densities['ho-165']
        er162 = densities['er-162']
        #er163 = densities['er-163']
        er164 = densities['er-164']
        #er165 = densities['er-165']
        er166 = densities['er-166']
        er167 = densities['er-167']
        er168 = densities['er-168']
        #er169 = densities['er-169']
        er170 = densities['er-170']
        lu175 = densities['lu-175']
        lu176 = densities['lu-176']
        hf174 = densities['hf-174']
        #hf175 = densities['hf-175']
        hf176 = densities['hf-176']
        hf177 = densities['hf-177']
        hf178 = densities['hf-178']
        hf179 = densities['hf-179']
        hf180 = densities['hf-180']
        ta181 = densities['ta-181']
        ta182 = densities['ta-182']
        w182 = densities['w-182']
        w183 = densities['w-183']
        w184 = densities['w-184']
        #w185 = densities['w-185']
        w186 = densities['w-186']
        re185 = densities['re-185']
        #re186 = densities['re-186']
        re187 = densities['re-187']
        ir191 = densities['ir-191']
        #ir192 = densities['ir-192']
        ir193 = densities['ir-193']
        #au193 = densities['au-193']
        #au194 = densities['au-194']
        #au195 = densities['au-195']
        #au196 = densities['au-196']
        au197 = densities['au-197']
        #au198 = densities['au-198']
        #au199 = densities['au-199']
        #au200 = densities['au-200']
        hg196 = densities['hg-196']
        #hg197 = densities['hg-197']
        hg198 = densities['hg-198']
        hg199 = densities['hg-199']
        hg200 = densities['hg-200']
        hg201 = densities['hg-201']
        hg202 = densities['hg-202']
        #hg203 = densities['hg-203']
        hg204 = densities['hg-204']
        pb204 = densities['pb-204']
        #pb205 = densities['pb-205']
        pb206 = densities['pb-206']
        pb207 = densities['pb-207']
        pb208 = densities['pb-208']
        bi209 = densities['bi-209']
        ra223 = densities['ra-223']
        ra224 = densities['ra-224']
        ra225 = densities['ra-225']
        ra226 = densities['ra-226']
        #ac224 = densities['ac-224']
        ac225 = densities['ac-225']
        ac226 = densities['ac-226']
        ac227 = densities['ac-227']
        #ac228 = densities['ac-228']
        u231 = densities['u-231']
        u232 = densities['u-232']
        u233 = densities['u-233']
        u237 = densities['u-237']
        u239 = densities['u-239']
        #th226 = densities['th-226']
        th227 = densities['th-227']
        th228 = densities['th-228']
        th229 = densities['th-229']
        th230 = densities['th-230']
        th231 = densities['th-231']
        th232 = densities['th-232']
        th233 = densities['th-233']
        th234 = densities['th-234']
        #pa228 = densities['pa-228']
        pa229 = densities['pa-229']
        pa230 = densities['pa-230']
        pa231 = densities['pa-231']
        pa232 = densities['pa-232']
        pa233 = densities['pa-233']
        #pa234 = densities['pa-234']
        #pa235 = densities['pa-235']
        pu236 = densities['pu-236']
        pu237 = densities['pu-237']
        pu238 = densities['pu-238']
        pu239 = densities['pu-239']
        pu240 = densities['pu-240']
        pu241 = densities['pu-241']
        pu242 = densities['pu-242']
        pu243 = densities['pu-243']
        pu244 = densities['pu-244']
        #pu245 = densities['pu-245']
        pu246 = densities['pu-246']
        #pu247 = densities['pu-247']
        np234 = densities['np-234']
        np235 = densities['np-235']
        np236 = densities['np-236']
        np237 = densities['np-237']
        np238 = densities['np-238']
        np239 = densities['np-239']
        #np240 = densities['np-240']
        #np241 = densities['np-241']
        cm240 = densities['cm-240']
        cm241 = densities['cm-241']
        cm242 = densities['cm-242']
        cm243 = densities['cm-243']
        cm244 = densities['cm-244']
        cm245 = densities['cm-245']
        cm246 = densities['cm-246']
        cm247 = densities['cm-247']
        cm248 = densities['cm-248']
        cm249 = densities['cm-249']
        cm250 = densities['cm-250']
        #cm251 = densities['cm-251']
        es251 = densities['es-251']
        es252 = densities['es-252']
        es253 = densities['es-253']
        es254 = densities['es-254']
        es255 = densities['es-255']
        #am239 = densities['am-239']
        am240 = densities['am-240']
        am241 = densities['am-241']
        am242 = densities['am-242']
        am243 = densities['am-243']
        am244 = densities['am-244']
        #am245 = densities['am-245']
        #am246 = densities['am-246']
        #am247 = densities['am-247']
        bk245 = densities['bk-245']
        bk246 = densities['bk-246']
        bk247 = densities['bk-247']
        bk248 = densities['bk-248']
        bk249 = densities['bk-249']
        bk250 = densities['bk-250']
        #bk251 = densities['bk-251']
        cf246 = densities['cf-246']
        cf248 = densities['cf-248']
        cf249 = densities['cf-249']
        cf250 = densities['cf-250']
        cf251 = densities['cf-251']
        cf252 = densities['cf-252']
        cf253 = densities['cf-253']
        cf254 = densities['cf-254']
        #cf255 = densities['cf-255']


        burned_salt_vector = [li6, li7, be9, f19, u234, u235, u236, u238, h1, h2, h3, he3, he4, be7, b10, b11, n14,
                            n15, o16, o17, na23, mg24, mg25, mg26, al27, si28, si29, si30, p31, s32, s33,
                            cl35, cl37, ar36, ar38, ar40, as74, as75, k39, k40,
                            k41, ca40, ca42, ca43, ca44, ca46, ca48, sc45, ti46, ti47, ti48, ti49,
                            ti50, cr50, cr52, cr53, cr54, mn55, fe54, fe56, fe57, fe58, co58, co59, ni58, ni59, ni60, ni61,
                            ni62, ni64, cu63, cu65, ga69, ga71, ge70, ge72, ge73, ge74, ge76, se74,
                            se76, se77, se78, se79, se80, se82, br79, br81, kr78, kr80, kr82, kr83,
                            kr84, kr85, kr86, rb85, rb86, rb87, sr84, sr86, sr87, sr88, sr89, sr90, y89, y90, y91,
                            zr90, zr91, zr92, zr93, zr94, zr95, zr96, nb93, nb94, nb95, mo92, mo94, mo95, mo96,
                            mo97, mo98, mo99, mo100, tc99, ru96, ru98, ru99, ru100, ru101, ru102, ru103, ru104, ru105,
                            ru106, rh103, rh105, pd102, pd104, pd105, pd106, pd107, pd108, pd110,
                            ag107, ag109, ag111, cd106, cd108, cd110, cd111, cd112,
                            cd113, cd114, cd116, in113, in115, sn112, sn113, sn114,
                            sn115, sn116, sn117, sn118, sn119, sn120, sn122, sn123, sn124, sn125, sn126, sb121,
                            sb123, sb124, sb125, sb126, te120, te122, te123, te124, te125, te126,
                            te128, te130, te132, i127, i129, i130, i131, i135,
                            xe123, xe124, xe126, xe128, xe129, xe130, xe131, xe132, xe133, xe134, xe135,
                            xe136, ce136, ce138, ce139, ce140, ce141, ce142, ce143, ce144, pr141,
                            pr142, pr143, ba130, ba132, ba133, ba134, ba135, ba136, ba137, ba138,
                            ba140, la138, la139, la140, nd142, nd143, nd144, nd145, nd146, nd147, sm149, sm150,
                            sm151, sm152, sm153, sm154, eu151, eu152, eu153, eu154, eu155, eu156, eu157, gd152, gd153,
                            gd154, gd155, gd156, gd157, gd158, gd160, tb159, tb160, dy156, dy158, 
                            dy160, dy161, dy162, dy163, dy164, ho165, er162, er164, er166, er167, er168,
                            er170, lu175, lu176, hf174, hf176, hf177, hf178, hf179, hf180, ta181, ta182,
                            w182, w183, w184, w186, re185, re187, ir191, ir193, au197, hg196, hg198, hg199, hg200, hg201, hg202,
                            hg204, pb204, pb206, pb207, pb208, bi209, ra223, ra224, ra225, ra226, ac225,
                            ac226, ac227, u231, u232, u233, u237, u239, th227, th228, th229, th230, th231,
                            th232, th233, th234, pa229, pa230, pa231, pa232, pa233, pu236, pu237,
                            pu238, pu239, pu240, pu241, pu242, pu243, pu244, pu246, np234, np235, np236,
                            np237, np238, np239, cm240, cm241, cm242, cm243, cm244, cm245, cm246, cm247,
                            cm248, cm249, cm250, es251, es252, es253, es254, es255, am240, am241, am242,
                            am243, am244, bk245, bk246, bk247, bk248, bk249, bk250, cf246,
                            cf248, cf249, cf250, cf251, cf252, cf253, cf254]   # atoms/barn-cm density vector
        
        atoms_cm3_vector = [i * 1e24 for i in burned_salt_vector]   # Converts atoms/barn-cm to atoms/cm3
        
        last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter-1))

        # Multiply by salt volume in cm3 to get number of atoms:

        if iter == 1:
            salt_volume = self.V0
            bs_atoms = [i * salt_volume for i in atoms_cm3_vector]
        else:
            os.chdir(last_dep_step_path)
            filename = "Salt_volume_dep_step_{}.out".format(iter-1)
            last_dep_step_vol = genfromtxt("{}".format(filename), delimiter='')   # MTiHM from the previous depletion step
            bs_atoms = [i * last_dep_step_vol for i in atoms_cm3_vector]  # Actual burned salt gram quantities after multiplying by previous depletion step MTiHM value
        
        return bs_atoms
    

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

    def get_refuel_atoms(self, rvol):
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
        return iso_atoms
    
    #####################################
    #    Mix Burned Salt w/ Refuel
    #####################################

    def mix_salts(self, rvol, iter):
        'Mix the burned salt from the previous step with the refuel salt.'

        # ***** Sum Burned Salt and Refuel Salt Atoms *****

        bs_atoms = self.get_burned_salt_atoms(iter)       # Number of atoms for each isotope in burned salt mat
        refuel_atoms = self.get_refuel_atoms(rvol)    # Number of atoms for each isotope in refuel salt mat

        # Burned Salt FLiBe Constituents:

        Li6_bs = bs_atoms[0]
        Li7_bs = bs_atoms[1]
        Be9_bs = bs_atoms[2]
        F19_bs = bs_atoms[3]
        U234_bs = bs_atoms[4]
        U235_bs = bs_atoms[5]
        U236_bs = bs_atoms[6]
        U238_bs = bs_atoms[7]
        h1 = bs_atoms[8]
        h2 = bs_atoms[9]
        h3 = bs_atoms[10]
        he3 = bs_atoms[11]
        he4 = bs_atoms[12]
        be7 = bs_atoms[13]
      #  be8 = bs_atoms[14]
        b10 = bs_atoms[14]
        b11 = bs_atoms[15]
        #b12 = bs_atoms[17]
        n14 = bs_atoms[16]
        n15 = bs_atoms[17]
        #n16 = bs_atoms[20]
        o16 = bs_atoms[18]
        o17 = bs_atoms[19]
        #f20 = bs_atoms[23]
        na23 = bs_atoms[20]
        mg24 = bs_atoms[21]
        mg25 = bs_atoms[22]
        mg26 = bs_atoms[23]
        al27 = bs_atoms[24]
        si28 = bs_atoms[25]
        si29 = bs_atoms[26]
        si30 = bs_atoms[27]
        p31 = bs_atoms[28]
        s32 = bs_atoms[29]
        s33 = bs_atoms[30]
       # s35 = bs_atoms[32]
        #s34 = bs_atoms[36]
        cl35 = bs_atoms[31]
        #cl36 = bs_atoms[38]
        cl37 = bs_atoms[32]
        #cl38 = bs_atoms[40]
        ar36 = bs_atoms[33]
        #ar37 = bs_atoms[42]
        ar38 = bs_atoms[34]
        #ar39 = bs_atoms[44]
        ar40 = bs_atoms[35]
        #ar41 = bs_atoms[46]
        #as72 = bs_atoms[47]
        #as73 = bs_atoms[48]
        as74 = bs_atoms[36]
        as75 = bs_atoms[37]
        #as76 = bs_atoms[51]
        #as77 = bs_atoms[52]
        #as78 = bs_atoms[53]
        #as79 = bs_atoms[54]
        #as80 = bs_atoms[55]
        #as81 = bs_atoms[56]
        #as82 = bs_atoms[57]
        #as83 = bs_atoms[58]
        #as84 = bs_atoms[59]
        #as85 = bs_atoms[60]
        #as86 = bs_atoms[61]
        #as87 = bs_atoms[62]
        #as88 = bs_atoms[63]
        #as89 = bs_atoms[64]
        #as90 = bs_atoms[65]
        k39 = bs_atoms[38]
        k40 = bs_atoms[39]
        k41 = bs_atoms[40]
        ca40 = bs_atoms[41]
        #ca41 = bs_atoms[70]
        ca42 = bs_atoms[42]
        ca43 = bs_atoms[43]
        ca44 = bs_atoms[44]
        #ca45 = bs_atoms[74]
        ca46 = bs_atoms[45]
        #ca47 = bs_atoms[76]
        ca48 = bs_atoms[46]
        #ca49 = bs_atoms[78]
        sc45 = bs_atoms[47]
        ti46 = bs_atoms[48]
        ti47 = bs_atoms[49]
        ti48 = bs_atoms[50]
        ti49 = bs_atoms[51]
        ti50 = bs_atoms[52]
        cr50 = bs_atoms[53]
        cr52 = bs_atoms[54]
        cr53 = bs_atoms[55]
        cr54 = bs_atoms[56]
        mn55 = bs_atoms[57]
        fe54 = bs_atoms[58]
        fe56 = bs_atoms[59]
        fe57 = bs_atoms[60]
        fe58 = bs_atoms[61]
        co58 = bs_atoms[62]
        co59 = bs_atoms[63]
        ni58 = bs_atoms[64]
        ni59 = bs_atoms[65]
        ni60 = bs_atoms[66]
        ni61 = bs_atoms[67]
        ni62 = bs_atoms[68]
        ni64 = bs_atoms[69]
        cu63 = bs_atoms[70]
        #cu64 = bs_atoms[103]
        cu65 = bs_atoms[71]
        ga69 = bs_atoms[72]
        ga71 = bs_atoms[73]
        ge70 = bs_atoms[74]
        #ge71 = bs_atoms[108]
        ge72 = bs_atoms[75]
        ge73 = bs_atoms[76]
        ge74 = bs_atoms[77]
        #ge75 = bs_atoms[112]
        ge76 = bs_atoms[78]
        se74 = bs_atoms[79]
        #se75 = bs_atoms[115]
        se76 = bs_atoms[80]
        se77 = bs_atoms[81]
        se78 = bs_atoms[82]
        se79 = bs_atoms[83]
        se80 = bs_atoms[84]
        #se81 = bs_atoms[121]
        se82 = bs_atoms[85]
        br79 = bs_atoms[86]
        #br80 = bs_atoms[124]
        br81 = bs_atoms[87]
        kr78 = bs_atoms[88]
        #kr79 = bs_atoms[127]
        kr80 = bs_atoms[89]
        #kr81 = bs_atoms[129]
        kr82 = bs_atoms[90]
        kr83 = bs_atoms[91]
        kr84 = bs_atoms[92]
        kr85 = bs_atoms[93]
        kr86 = bs_atoms[94]
        rb85 = bs_atoms[95]
        rb86 = bs_atoms[96]
        rb87 = bs_atoms[97]
        sr84 = bs_atoms[98]
        #sr85 = bs_atoms[139]
        sr86 = bs_atoms[99]
        sr87 = bs_atoms[100]
        sr88 = bs_atoms[101]
        sr89 = bs_atoms[102]
        sr90 = bs_atoms[103]
        y89 = bs_atoms[104]
        y90 = bs_atoms[105]
        y91 = bs_atoms[106]
        #zr89 = bs_atoms[148]
        zr90 = bs_atoms[107]
        zr91 = bs_atoms[108]
        zr92 = bs_atoms[109]
        zr93 = bs_atoms[110]
        zr94 = bs_atoms[111]
        zr95 = bs_atoms[112]
        zr96 = bs_atoms[113]
        nb93 = bs_atoms[114]
        nb94 = bs_atoms[115]
        nb95 = bs_atoms[116]
        mo92 = bs_atoms[117]
        #mo93 = bs_atoms[160]
        mo94 = bs_atoms[118]
        mo95 = bs_atoms[119]
        mo96 = bs_atoms[120]
        mo97 = bs_atoms[121]
        mo98 = bs_atoms[122]
        mo99 = bs_atoms[123]
        mo100 = bs_atoms[124]
        tc99 = bs_atoms[125]
        ru96 = bs_atoms[126]
        #ru97 = bs_atoms[170]
        ru98 = bs_atoms[127]
        ru99 = bs_atoms[128]
        ru100 = bs_atoms[129]
        ru101 = bs_atoms[130]
        ru102 = bs_atoms[131]
        ru103 = bs_atoms[132]
        ru104 = bs_atoms[133]
        ru105 = bs_atoms[134]
        ru106 = bs_atoms[135]
        rh103 = bs_atoms[136]
        #rh104 = bs_atoms[181]
        rh105 = bs_atoms[137]
        pd102 = bs_atoms[138]
        #pd103 = bs_atoms[184]
        pd104 = bs_atoms[139]
        pd105 = bs_atoms[140]
        pd106 = bs_atoms[141]
        pd107 = bs_atoms[142]
        pd108 = bs_atoms[143]
        #pd109 = bs_atoms[190]
        pd110 = bs_atoms[144]
        ag107 = bs_atoms[145]
        #ag108 = bs_atoms[193]
        ag109 = bs_atoms[146]
        #ag110 = bs_atoms[195]
        ag111 = bs_atoms[147]
        #cd105 = bs_atoms[197]
        cd106 = bs_atoms[148]
        #cd107 = bs_atoms[199]
        cd108 = bs_atoms[149]
        #cd109 = bs_atoms[201]
        cd110 = bs_atoms[150]
        cd111 = bs_atoms[151]
        cd112 = bs_atoms[152]
        cd113 = bs_atoms[153]
        cd114 = bs_atoms[154]
        #cd115 = bs_atoms[207]
        cd116 = bs_atoms[155]
        #cd117 = bs_atoms[209]
        #cd118 = bs_atoms[210]
        #cd119 = bs_atoms[211]
        #cd120 = bs_atoms[212]
        #cd121 = bs_atoms[213]
        #cd122 = bs_atoms[214]
        #cd123 = bs_atoms[215]
        #cd124 = bs_atoms[216]
        #cd125 = bs_atoms[217]
        #cd126 = bs_atoms[218]
        #cd127 = bs_atoms[219]
        #cd128 = bs_atoms[220]
        #cd129 = bs_atoms[221]
        #cd130 = bs_atoms[222]
        #cd131 = bs_atoms[223]
        #cd132 = bs_atoms[224]
        in113 = bs_atoms[156]
        #in114 = bs_atoms[226]
        in115 = bs_atoms[157]
        sn112 = bs_atoms[158]
        sn113 = bs_atoms[159]
        sn114 = bs_atoms[160]
        sn115 = bs_atoms[161]
        sn116 = bs_atoms[162]
        sn117 = bs_atoms[163]
        sn118 = bs_atoms[164]
        sn119 = bs_atoms[165]
        sn120 = bs_atoms[166]
        #sn121 = bs_atoms[237]
        sn122 = bs_atoms[167]
        sn123 = bs_atoms[168]
        sn124 = bs_atoms[169]
        sn125 = bs_atoms[170]
        sn126 = bs_atoms[171]
        sb121 = bs_atoms[172]
        #sb122 = bs_atoms[244]
        sb123 = bs_atoms[173]
        sb124 = bs_atoms[174]
        sb125 = bs_atoms[175]
        sb126 = bs_atoms[176]
        te120 = bs_atoms[177]
        #te121 = bs_atoms[250]
        te122 = bs_atoms[178]
        te123 = bs_atoms[179]
        te124 = bs_atoms[180]
        te125 = bs_atoms[181]
        te126 = bs_atoms[182]
        #te127 = bs_atoms[256]
        te128 = bs_atoms[183]
        #te129 = bs_atoms[258]
        te130 = bs_atoms[184]
        #te131 = bs_atoms[260]
        te132 = bs_atoms[185]
        i127 = bs_atoms[186]
        #i128 = bs_atoms[263]
        i129 = bs_atoms[187]
        i130 = bs_atoms[188]
        i131 = bs_atoms[189]
        #i132 = bs_atoms[267]
        #i133 = bs_atoms[268]
        #i134 = bs_atoms[269]
        i135 = bs_atoms[190]
        xe123 = bs_atoms[191]
        xe124 = bs_atoms[192]
        #xe125 = bs_atoms[273]
        xe126 = bs_atoms[193]
        #xe127 = bs_atoms[275]
        xe128 = bs_atoms[194]
        xe129 = bs_atoms[195]
        xe130 = bs_atoms[196]
        xe131 = bs_atoms[197]
        xe132 = bs_atoms[198]
        xe133 = bs_atoms[199]
        xe134 = bs_atoms[200]
        xe135 = bs_atoms[201]
        xe136 = bs_atoms[202]
        #ba128 = bs_atoms[285]
        #ba129 = bs_atoms[286]
        ba130 = bs_atoms[203]
        #ba131 = bs_atoms[288]
        ba132 = bs_atoms[204]
        ba133 = bs_atoms[205]
        ba134 = bs_atoms[206]
        ba135 = bs_atoms[207]
        ba136 = bs_atoms[208]
        ba137 = bs_atoms[209]
        ba138 = bs_atoms[210]
        #ba139 = bs_atoms[296]
        ba140 = bs_atoms[211]
        #ba141 = bs_atoms[298]
        #ba142 = bs_atoms[299]
        #ba143 = bs_atoms[300]
        #ba144 = bs_atoms[301]
        #ba145 = bs_atoms[302]
        #ba146 = bs_atoms[303]
        #ba147 = bs_atoms[304]
        #ba148 = bs_atoms[305]
        #ba149 = bs_atoms[306]
        #ba150 = bs_atoms[307]
        #ba151 = bs_atoms[308]
        #ba152 = bs_atoms[309]
        #ba153 = bs_atoms[310]
        la138 = bs_atoms[212]
        la139 = bs_atoms[213]
        la140 = bs_atoms[214]
        #ce134 = bs_atoms[314]
        #ce135 = bs_atoms[315]
        ce136 = bs_atoms[215]
        #ce137 = bs_atoms[218]
        ce138 = bs_atoms[216]
        ce139 = bs_atoms[217]
        ce140 = bs_atoms[218]
        ce141 = bs_atoms[219]
        ce142 = bs_atoms[220]
        ce143 = bs_atoms[221]
        ce144 = bs_atoms[222]
        #ce145 = bs_atoms[325]
        #ce146 = bs_atoms[326]
        #ce147 = bs_atoms[327]
        #ce148 = bs_atoms[328]
        #ce149 = bs_atoms[329]
        #ce150 = bs_atoms[330]
        #ce151 = bs_atoms[331]
        #ce152 = bs_atoms[332]
        #ce153 = bs_atoms[333]
        #ce154 = bs_atoms[334]
        #ce155 = bs_atoms[335]
        #ce156 = bs_atoms[336]
        #ce157 = bs_atoms[337]
        pr141 = bs_atoms[223]
        pr142 = bs_atoms[224]
        pr143 = bs_atoms[225]
        nd142 = bs_atoms[226]
        nd143 = bs_atoms[227]
        nd144 = bs_atoms[228]
        nd145 = bs_atoms[229]
        nd146 = bs_atoms[230]
        nd147 = bs_atoms[231]
        sm149 = bs_atoms[232]
        sm150 = bs_atoms[233]
        sm151 = bs_atoms[234]
        sm152 = bs_atoms[235]
        sm153 = bs_atoms[236]
        sm154 = bs_atoms[237]
        eu151 = bs_atoms[238]
        eu152 = bs_atoms[239]
        eu153 = bs_atoms[240]
        eu154 = bs_atoms[241]
        eu155 = bs_atoms[242]
        eu156 = bs_atoms[243]
        eu157 = bs_atoms[244]
        gd152 = bs_atoms[245]
        gd153 = bs_atoms[246]
        gd154 = bs_atoms[247]
        gd155 = bs_atoms[248]
        gd156 = bs_atoms[249]
        gd157 = bs_atoms[250]
        gd158 = bs_atoms[251]
        #gd159 = bs_atoms[367]
        gd160 = bs_atoms[252]
        tb159 = bs_atoms[253]
        tb160 = bs_atoms[254]
        dy156 = bs_atoms[255]
        #dy157 = bs_atoms[372]
        dy158 = bs_atoms[256]
        #dy159 = bs_atoms[374]
        dy160 = bs_atoms[257]
        dy161 = bs_atoms[258]
        dy162 = bs_atoms[259]
        dy163 = bs_atoms[260]
        dy164 = bs_atoms[261]
        ho165 = bs_atoms[262]
        er162 = bs_atoms[263]
        #er163 = bs_atoms[382]
        er164 = bs_atoms[264]
        #er165 = bs_atoms[384]
        er166 = bs_atoms[265]
        er167 = bs_atoms[266]
        er168 = bs_atoms[267]
        #er169 = bs_atoms[388]
        er170 = bs_atoms[268]
        lu175 = bs_atoms[269]
        lu176 = bs_atoms[270]
        hf174 = bs_atoms[271]
        #hf175 = bs_atoms[393]
        hf176 = bs_atoms[272]
        hf177 = bs_atoms[273]
        hf178 = bs_atoms[274]
        hf179 = bs_atoms[275]
        hf180 = bs_atoms[276]
        ta181 = bs_atoms[277]
        ta182 = bs_atoms[278]
        w182 = bs_atoms[279]
        w183 = bs_atoms[280]
        w184 = bs_atoms[281]
        #w185 = bs_atoms[404]
        w186 = bs_atoms[282]
        re185 = bs_atoms[283]
        #re186 = bs_atoms[407]
        re187 = bs_atoms[284]
        ir191 = bs_atoms[285]
        #ir192 = bs_atoms[410]
        ir193 = bs_atoms[286]
        #au193 = bs_atoms[412]
        #au194 = bs_atoms[413]
        #au195 = bs_atoms[414]
        #au196 = bs_atoms[415]
        au197 = bs_atoms[287]
        #au198 = bs_atoms[417]
        #au199 = bs_atoms[418]
        #au200 = bs_atoms[419]
        hg196 = bs_atoms[288]
        #hg197 = bs_atoms[421]
        hg198 = bs_atoms[289]
        hg199 = bs_atoms[290]
        hg200 = bs_atoms[291]
        hg201 = bs_atoms[292]
        hg202 = bs_atoms[293]
        #hg203 = bs_atoms[427]
        hg204 = bs_atoms[294]
        pb204 = bs_atoms[295]
        #pb205 = bs_atoms[430]
        pb206 = bs_atoms[296]
        pb207 = bs_atoms[297]
        pb208 = bs_atoms[298]
        bi209 = bs_atoms[299]
        ra223 = bs_atoms[300]
        ra224 = bs_atoms[301]
        ra225 = bs_atoms[302]
        ra226 = bs_atoms[303]
        #ac224 = bs_atoms[439]
        ac225 = bs_atoms[304]
        ac226 = bs_atoms[305]
        ac227 = bs_atoms[306]
        #ac228 = bs_atoms[443]
        u231 = bs_atoms[307]
        u232 = bs_atoms[308]
        u233 = bs_atoms[309]
        u237 = bs_atoms[310]
        u239 = bs_atoms[311]
        #th226 = bs_atoms[449]
        th227 = bs_atoms[312]
        th228 = bs_atoms[313]
        th229 = bs_atoms[314]
        th230 = bs_atoms[315]
        th231 = bs_atoms[316]
        th232 = bs_atoms[317]
        th233 = bs_atoms[318]
        th234 = bs_atoms[319]
        #pa228 = bs_atoms[458]
        pa229 = bs_atoms[320]
        pa230 = bs_atoms[321]
        pa231 = bs_atoms[322]
        pa232 = bs_atoms[323]
        pa233 = bs_atoms[324]
        #pa234 = bs_atoms[464]
        #pa235 = bs_atoms[465]
        pu236 = bs_atoms[325]
        pu237 = bs_atoms[326]
        pu238 = bs_atoms[327]
        pu239 = bs_atoms[328]
        pu240 = bs_atoms[329]
        pu241 = bs_atoms[330]
        pu242 = bs_atoms[331]
        pu243 = bs_atoms[332]
        pu244 = bs_atoms[333]
        #pu245 = bs_atoms[474]
        pu246 = bs_atoms[334]
        #pu247 = bs_atoms[476]
        #pu238 = bs_atoms[336]
        np234 = bs_atoms[335]
        np235 = bs_atoms[336]
        np236 = bs_atoms[337]
        np237 = bs_atoms[338]
        np238 = bs_atoms[339]
        np239 = bs_atoms[340]
        #np240 = bs_atoms[484]
        #np241 = bs_atoms[485]
        cm240 = bs_atoms[341]
        cm241 = bs_atoms[342]
        cm242 = bs_atoms[343]
        cm243 = bs_atoms[344]
        cm244 = bs_atoms[345]
        cm245 = bs_atoms[346]
        cm246 = bs_atoms[347]
        cm247 = bs_atoms[348]
        cm248 = bs_atoms[349]
        cm249 = bs_atoms[350]
        cm250 = bs_atoms[351]
        #cm251 = bs_atoms[497]
        es251 = bs_atoms[352]
        es252 = bs_atoms[353]
        es253 = bs_atoms[354]
        es254 = bs_atoms[355]
        es255 = bs_atoms[356]
        #am239 = bs_atoms[503]
        am240 = bs_atoms[357]
        am241 = bs_atoms[358]
        am242 = bs_atoms[359]
        am243 = bs_atoms[360]
        am244 = bs_atoms[361]
        #am245 = bs_atoms[509]
        #am246 = bs_atoms[510]
        #am247 = bs_atoms[511]
        bk245 = bs_atoms[362]
        bk246 = bs_atoms[363]
        bk247 = bs_atoms[364]
        bk248 = bs_atoms[365]
        bk249 = bs_atoms[366]
        bk250 = bs_atoms[367]
        #bk251 = bs_atoms[518]
        cf246 = bs_atoms[368]
        cf248 = bs_atoms[369]
        cf249 = bs_atoms[370]
        cf250 = bs_atoms[371]
        cf251 = bs_atoms[372]
        cf252 = bs_atoms[373]
        cf253 = bs_atoms[374]
        cf254 = bs_atoms[375]
        #cf255 = bs_atoms[527]


        # Refuel Salt FLiBe SConstituents:

        Li6_rf = refuel_atoms[0]
        Li7_rf = refuel_atoms[1]
        Be9_rf = refuel_atoms[2]
        F19_rf = refuel_atoms[3]
        U234_rf = refuel_atoms[4]
        U235_rf = refuel_atoms[5]
        U236_rf = refuel_atoms[6]
        U238_rf = refuel_atoms[7]

        # Totals:

        Li6_tot = Li6_bs + Li6_rf
        Li7_tot = Li7_bs + Li7_rf
        Be9_tot = Be9_bs + Be9_rf
        F19_tot = F19_bs + F19_rf
        U234_tot = U234_bs + U234_rf
        U235_tot = U235_bs + U235_rf
        U236_tot = U236_bs + U236_rf
        U238_tot = U238_bs + U238_rf

        # ***** New Total Fuel Salt Mixed Atom Matrix *****

        mixed_atom_totals = [Li6_tot, Li7_tot, Be9_tot, F19_tot, U234_tot, U235_tot, U236_tot, U238_tot, h1, h2, h3, he3,
                             he4, be7, b10, b11, n14, n15, o16, o17, na23, mg24, mg25, mg26, al27, si28,
                             si29, si30, p31, s32, s33, cl35, cl37, ar36, ar38, ar40,
                             as74, as75, k39, k40, k41, ca40, ca42, ca43, ca44, ca46, ca48, sc45,
                             ti46, ti47, ti48, ti49, ti50, cr50, cr52, cr53, cr54, mn55, fe54, fe56, fe57, fe58, co58, co59,
                             ni58, ni59, ni60, ni61, ni62, ni64, cu63, cu65, ga69, ga71, ge70, ge72, ge73, ge74,
                             ge76, se74, se76, se77, se78, se79, se80, se82, br79, br81, kr78,
                             kr80, kr82, kr83, kr84, kr85, kr86, rb85, rb86, rb87, sr84, sr86, sr87, sr88, sr89,
                             sr90, y89, y90, y91, zr90, zr91, zr92, zr93, zr94, zr95, zr96, nb93, nb94, nb95, mo92,
                             mo94, mo95, mo96, mo97, mo98, mo99, mo100, tc99, ru96, ru98, ru99, ru100, ru101, ru102,
                             ru103, ru104, ru105, ru106, rh103, rh105, pd102, pd104, pd105, pd106, pd107, pd108,
                             pd110, ag107, ag109, ag111, cd106, cd108, cd110, cd111, cd112,
                             cd113, cd114, cd116, in113, in115, sn112, sn113, sn114,
                             sn115, sn116, sn117, sn118, sn119, sn120, sn122, sn123, sn124, sn125, sn126, sb121,
                             sb123, sb124, sb125, sb126, te120, te122, te123, te124, te125, te126,
                             te128, te130, te132, i127, i129, i130, i131, i135,
                             xe123, xe124, xe126, xe128, xe129, xe130, xe131, xe132, xe133, xe134, xe135,
                             xe136, ce136, ce138, ce139, ce140, ce141, ce142, ce143, ce144, pr141,
                             pr142, pr143, ba130, ba132, ba133, ba134, ba135, ba136, ba137, ba138,
                             ba140, la138, la139, la140, nd142, nd143, nd144, nd145, nd146, nd147, sm149, sm150,
                             sm151, sm152, sm153, sm154, eu151, eu152, eu153, eu154, eu155, eu156, eu157, gd152, gd153,
                             gd154, gd155, gd156, gd157, gd158, gd160, tb159, tb160, dy156, dy158,
                             dy160, dy161, dy162, dy163, dy164, ho165, er162, er164, er166, er167, er168,
                             er170, lu175, lu176, hf174, hf176, hf177, hf178, hf179, hf180, ta181, ta182,
                             w182, w183, w184, w186, re185, re187, ir191, ir193,
                             au197, hg196, hg198, hg199, hg200, hg201, hg202,
                             hg204, pb204, pb206, pb207, pb208, bi209, ra223, ra224, ra225, ra226, ac225,
                             ac226, ac227, u231, u232, u233, u237, u239, th227, th228, th229, th230, th231,
                             th232, th233, th234, pa229, pa230, pa231, pa232, pa233, pu236, pu237,
                             pu238, pu239, pu240, pu241, pu242, pu243, pu244, pu246, np234, np235, np236,
                             np237, np238, np239, cm240, cm241, cm242, cm243, cm244, cm245, cm246, cm247,
                             cm248, cm249, cm250, es251, es252, es253, es254, es255, am240, am241, am242,
                             am243, am244, bk245, bk246, bk247, bk248, bk249, bk250, cf246,
                             cf248, cf249, cf250, cf251, cf252, cf253, cf254] 

        last_dep_step = iter - 1
        last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter-1))
        filename = "Salt_volume_dep_step_{}.out".format(last_dep_step)
        if iter == 1:
            bs_vol = self.V0     # Initial volume of salt in cm3
        else:
            os.chdir(last_dep_step_path)
            bs_vol = genfromtxt("{}".format(filename), delimiter='')   # Total fuel salt volume for previous depletion step

        tot_salt_vol = bs_vol + rvol

        # Finding atom densities of each from total salt volume (tot_salt_vol):
        mixed_atom_den_cm3 = [i / tot_salt_vol for i in mixed_atom_totals]    # Mixed salt atom density in atoms/cm3 for each isotope
        mixed_atom_den = [i * 1e-24 for i in mixed_atom_den_cm3]       # Mixed salt atom density vector in atoms/barn-cm for each isotope
        return mixed_atom_den
    
    def write_scale_mat(self, rvol, iter):
        'Returns a SCALE material composition block for the refuel salt mixed in with the burned salt.'
        enrich_percent = self.renrich*100    # Enrichment percent for refuel
        mixed_salt_aden = self.mix_salts(rvol, iter)
        isotopes = ["li-6", "li-7", "be-9", "f-19", "u-234", "u-235", "u-236", "u-238", "h-1", "h-2", "h-3", "he-3", "he-4", "be-7",
                    "b-10", "b-11", "n-14", "n-15", "o-16", "o-17", "na-23", "mg-24", "mg-25", "mg-26",
                    "al-27", "si-28", "si-29", "si-30", "p-31", "s-32", "s-33", "cl-35", "cl-37",
                    "ar-36", "ar-38", "ar-40", "as-74", "as-75", "k-39",
                    "k-40", "k-41", "ca-40", "ca-42", "ca-43", "ca-44", "ca-46", "ca-48", "sc-45",
                    "ti-46", "ti-47", "ti-48", "ti-49", "ti-50", "cr-50", "cr-52", "cr-53", "cr-54", "mn-55", "fe-54", "fe-56", "fe-57", "fe-58", 
                    "co-58", "co-59", "ni-58", "ni-59", "ni-60", "ni-61", "ni-62", "ni-64", "cu-63", "cu-65", "ga-69", "ga-71", "ge-70", 
                    "ge-72", "ge-73", "ge-74", "ge-76", "se-74", "se-76", "se-77", "se-78", "se-79", "se-80",
                    "se-82", "br-79", "br-81", "kr-78", "kr-80", "kr-82", "kr-83", "kr-84", "kr-85", "kr-86", "rb-85",
                    "rb-86", "rb-87", "sr-84", "sr-86", "sr-87", "sr-88", "sr-89", "sr-90", "y-89", "y-90", "y-91", "zr-90", "zr-91",
                    "zr-92", "zr-93", "zr-94", "zr-95", "zr-96", "nb-93", "nb-94", "nb-95", "mo-92", "mo-94", "mo-95", "mo-96", "mo-97", "mo-98",
                    "mo-99", "mo-100", "tc-99", "ru-96", "ru-98", "ru-99", "ru-100", "ru-101", "ru-102", "ru-103", "ru-104", "ru-105", "ru-106",
                    "rh-103", "rh-105", "pd-102", "pd-104", "pd-105", "pd-106", "pd-107", "pd-108", "pd-110", "ag-107",
                    "ag-109", "ag-111", "cd-106", "cd-108", "cd-110", "cd-111", "cd-112", "cd-113",
                    "cd-114", "cd-116", "in-113", "in-115", "sn-112", "sn-113", "sn-114", "sn-115",
                    "sn-116", "sn-117", "sn-118", "sn-119", "sn-120", "sn-122", "sn-123", "sn-124", "sn-125", "sn-126", "sb-121",
                    "sb-123", "sb-124", "sb-125", "sb-126", "te-120", "te-122", "te-123", "te-124", "te-125", "te-126", "te-128",
                    "te-130", "te-132", "i-127", "i-129", "i-130", "i-131", "i-135", "xe-123",
                    "xe-124", "xe-126", "xe-128", "xe-129", "xe-130", "xe-131", "xe-132", "xe-133", "xe-134", "xe-135", "xe-136",
                    "ce-136", "ce-138", "ce-139", "ce-140", "ce-141", "ce-142", "ce-143", "ce-144", "pr-141", "pr-142",
                    "pr-143", "ba-130", "ba-132", "ba-133", "ba-134", "ba-135", "ba-136", "ba-137", "ba-138",
                    "ba-140", "la-138", "la-139", "la-140", "nd-142", "nd-143", "nd-144", "nd-145", "nd-146", "nd-147", "sm-149", "sm-150", "sm-151",
                    "sm-152", "sm-153", "sm-154", "eu-151", "eu-152", "eu-153", "eu-154", "eu-155", "eu-156", "eu-157", "gd-152", "gd-153", "gd-154",
                    "gd-155", "gd-156", "gd-157", "gd-158", "gd-160", "tb-159", "tb-160", "dy-156", "dy-158", "dy-160",
                    "dy-161", "dy-162", "dy-163", "dy-164", "ho-165", "er-162", "er-164", "er-166", "er-167", "er-168",
                    "er-170", "lu-175", "lu-176", "hf-174", "hf-176", "hf-177", "hf-178", "hf-179", "hf-180", "ta-181", "ta-182", "w-182",
                    "w-183", "w-184", "w-186", "re-185", "re-187", "ir-191", "ir-193", "au-197", "hg-196", "hg-198", "hg-199", "hg-200", "hg-201", "hg-202",
                    "hg-204", "pb-204", "pb-206", "pb-207", "pb-208", "bi-209", "ra-223", "ra-224", "ra-225", "ra-226", "ac-225",
                    "ac-226", "ac-227", "u-231", "u-232", "u-233", "u-237", "u-239", "th-227", "th-228", "th-229", "th-230", "th-231",
                    "th-232", "th-233", "th-234", "pa-229", "pa-230", "pa-231", "pa-232", "pa-233", "pu-236", "pu-237",
                    "pu-238", "pu-239", "pu-240", "pu-241", "pu-242", "pu-243", "pu-244", "pu-246", "np-234", "np-235", "np-236",
                    "np-237", "np-238", "np-239", "cm-240", "cm-241", "cm-242", "cm-243", "cm-244", "cm-245", "cm-246", "cm-247",
                    "cm-248", "cm-249", "cm-250", "es-251", "es-252", "es-253", "es-254", "es-255", "am-240", "am-241", "am-242",
                    "am-243", "am-244", "bk-245", "bk-246", "bk-247", "bk-248", "bk-249", "bk-250", "cf-246",
                    "cf-248", "cf-249", "cf-250", "cf-251", "cf-252", "cf-253", "cf-254"]
        scale_mat = "' Burned EIRENE fuel salt " + "and mixed refuel with enrichment of " + str(enrich_percent) + "%" + "\n"
        for i in range(376):
            scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], mixed_salt_aden[i])
        return scale_mat

    
        ##################################
        #     Get New MTiHM Value
        ##################################

    def get_burned_salt_MTHM(self):
        'Returns the masses in grams of MTHM actinides (excluding Actinium - SCALE does not consider Actinium as part of MTHM) in the burned salt.'
        output = subprocess.run([f"{SCALE_bin_path}/obiwan", "view", "-format=csv", "-prec=10", "-units=gram", "-idform='{:Ee}{:AAA}{:m}'", self.f71_name], capture_output=True)
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
                densities[nuclide] = float(data[2])  # The [2] here is what causes the code to return only the densities at position 2 of the f71 file (position 0 is isotope names, 1 is BOC data)
        # **** Sum all Actinide HM masses: ****
        u231 = densities['u-231']
        u232 = densities['u-232']
        u233 = densities['u-233']
        u234 = densities['u-234']
        u235 = densities['u-235']
        u235m = densities['u-235m']
        u236 = densities['u-236']
        u237 = densities['u-237']
        u238 = densities['u-238']
        u239 = densities['u-239']
        th226 = densities['th-226']
        th227 = densities['th-227']
        th228 = densities['th-228']
        th229 = densities['th-229']
        th230 = densities['th-230']
        th231 = densities['th-231']
        th232 = densities['th-232']
        th233 = densities['th-233']
        th234 = densities['th-234']
        pa228 = densities['pa-228']
        pa229 = densities['pa-229']
        pa230 = densities['pa-230']
        pa231 = densities['pa-231']
        pa232 = densities['pa-232']
        pa233 = densities['pa-233']
        pa234 = densities['pa-234']
        pa234m = densities['pa-234m']
        pa235 = densities['pa-235']
        pu236 = densities['pu-236']
        pu237 = densities['pu-237']
        pu237m = densities['pu-237m']
        pu238 = densities['pu-238']
        pu239 = densities['pu-239']
        pu240 = densities['pu-240']
        pu241 = densities['pu-241']
        pu242 = densities['pu-242']
        pu243 = densities['pu-243']
        pu244 = densities['pu-244']
        pu245 = densities['pu-245']
        pu246 = densities['pu-246']
        pu247 = densities['pu-247']
        np234 = densities['np-234']
        np235 = densities['np-235']
        np236 = densities['np-236']
        np236m = densities['np-236m']
        np237 = densities['np-237']
        np238 = densities['np-238']
        np239 = densities['np-239']
        np240 = densities['np-240']
        np240m = densities['np-240m']
        np241 = densities['np-241']
        am239 = densities['am-239']
        am240 = densities['am-240']
        am241 = densities['am-241']
        am242 = densities['am-242']
        am242m = densities['am-242m']
        am243 = densities['am-243']
        am244 = densities['am-244']
        am244m = densities['am-244m']
        am245 = densities['am-245']
        am246 = densities['am-246']
        am246m = densities['am-246m']
        am247 = densities['am-247']
        cm240 = densities['cm-240']
        cm241 = densities['cm-241']
        cm242 = densities['cm-242']
        cm243 = densities['cm-243']
        cm244 = densities['cm-244']
        cm245 = densities['cm-245']
        cm246 = densities['cm-246']
        cm247 = densities['cm-247']
        cm248 = densities['cm-248']
        cm249 = densities['cm-249']
        cm250 = densities['cm-250']
        cm251 = densities['cm-251']
        es251 = densities['es-251']
        es252 = densities['es-252']
        es253 = densities['es-253']
        es254 = densities['es-254']
        es254m = densities['es-254m']
        es255 = densities['es-255']
        bk245 = densities['bk-245']
        bk246 = densities['bk-246']
        bk247 = densities['bk-247']
        bk248 = densities['bk-248']
        bk249 = densities['bk-249']
        bk250 = densities['bk-250']
        bk251 = densities['bk-251']
        cf246 = densities['cf-246']
        cf248 = densities['cf-248']
        cf249 = densities['cf-249']
        cf250 = densities['cf-250']
        cf251 = densities['cf-251']
        cf252 = densities['cf-252']
        cf253 = densities['cf-253']
        cf254 = densities['cf-254']
        cf255 = densities['cf-255']


        HMs = (u231 + u232 + u233 + u234 + u235 + u235m + u236 + u237 + u238 + u239 + th226 + th227 + th228 + th229 + th230 + th231 + th232
        + th233 + th234 + pa228 + pa229 + pa230 + pa231 + pa232 + pa233 + pa234 + pa234m + pa235 + pu236 + pu237 + pu237m + pu238 + pu239
        + pu240 + pu241 + pu242 + pu243 + pu244 + pu245 + pu246 + pu247 + np234 + np235 + np236 + np236m + np237 + np238 + np239 + np240
        + np240m + np241 + am239 + am240 + am241 + am242 + am242m + am243 + am244 + am244m + am245 + am246 + am246m + am247 + cm240 + cm241
        + cm242 + cm243 + cm244 + cm245 + cm246 + cm247 + cm248 + cm249 + cm250 + cm251 + es251 + es252 + es253 + es254 + es254m + es255
        + bk245 + bk246 + bk247 + bk248 + bk249 + bk250 + bk251 + cf246 + cf248 + cf249 + cf250 + cf251 + cf252 + cf253 + cf254 + cf255)   # Total heavy metal mass in grams

        MTHMs = HMs*1e-6     # Convert HMs mass in grams to MTHMs

        return MTHMs
    
    def get_refuel_MTHM(self, iter):
        refuel_wf = self.get_refuel_wf()             # Isotopic weight fraction list for refuel salt
        refuel_dens = self.get_refuel_den()          # Density of refuel salt in g/cm3
        rvol = self.get_crit_refuel(iter)            # Critical refuel amount added in cm3
        total_refuel_mass = rvol*refuel_dens         # Total mass of refuel salt in grams
        # Calculate mass of each isotope in the refuel:
        iso_mass = [i * total_refuel_mass for i in refuel_wf]
        u234 = iso_mass[4]    # U-234 mass in grams
        u235 = iso_mass[5]    # U-235 mass in grams
        u236 = iso_mass[6]    # U-236 mass in grams
        u238 = iso_mass[7]    # U-238 mass in grams
        refuel_HM = (u234 + u235 + u236 + u238)      # Total refuel mass of HMs in grams
        refuel_MTHM = refuel_HM*1e-6                 # Total refuel MTHM
        return refuel_MTHM
    
    def write_MTiHM_file(self,iter):
        'Finds the new MTiHM value for the current depletion step and saves it to a .out file to be read in the next depletion step.'
        # ADD: Something that references self.vol_doubled variable and checks whether it is True or False. If True, it halves the MTiHM value.
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        if iter == 1:
            os.chdir(self.BOC_f71_path)
        else:
            last_dep_step = iter - 1
            last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(last_dep_step))
            os.chdir(last_dep_step_path)
        bs_MTHM = self.get_burned_salt_MTHM()
        os.chdir(main_path)
        refuel_MTHM = self.get_refuel_MTHM(iter)

        if iter == 1:
            true_bs_MTHM = bs_MTHM*self.init_MTiHM
            os.chdir(main_path)
            tot_MTiHM = true_bs_MTHM + refuel_MTHM
            f = open("MTiHM_dep_step_{}.out".format(iter), 'w')
            f.write("{}".format(tot_MTiHM))
            f.close()
        else:
            os.chdir(last_dep_step_path)
            filename = "MTiHM_dep_step_{}.out".format(iter-1)
            data = genfromtxt("{}".format(filename), delimiter='')
            true_bs_MTHM = bs_MTHM*data
            tot_MTiHM = true_bs_MTHM + refuel_MTHM
            if self.vol_doubled == True:
                new_tot_MTiHM = tot_MTiHM/2
            else:
                new_tot_MTiHM = tot_MTiHM
            os.chdir(main_path)
            f = open("MTiHM_dep_step_{}.out".format(iter), 'w')
            f.write("{}".format(new_tot_MTiHM))
            f.close()

    #################################
    #        Add in Refuel
    #################################
            
    def write_fuel_salt_volume(self, rvol, iter, main_path):
        'Writes a .out file with the total fuel salt volume for the current depletion step.'
        current_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter-1))
        self.vol_doubled = False  # Check parameter; changed to True if the volume has doubled. Used for halving the MTiHM value if the volume has doubled (add line to MTiHM file writing function that does this.)
        if iter == 1:
            volume = self.V0
            total_vol = volume + rvol    # Total fuel salt volume in cm^3
            if total_vol >= 2*self.V0:
                self.vol_doubled = True
                new_vol = total_vol/2
            else:
                new_vol = total_vol
            f = open("Salt_volume_dep_step_{}.out".format(iter), 'w')
            f.write("{}".format(new_vol))
            f.close()
        else:
            os.chdir(last_dep_step_path)
            filename = "Salt_volume_dep_step_{}.out".format(iter-1)
            volume = genfromtxt("{}".format(filename), delimiter='')
            total_vol = volume + rvol
            if total_vol >= 2*self.V0:
                self.vol_doubled = True
                new_vol = total_vol/2
            else:
                new_vol = total_vol
              
            os.chdir(current_dep_step_path)
            os.chdir(main_path)
            f = open("Salt_volume_dep_step_{}.out".format(iter), 'w')
            f.write("{}".format(new_vol))
            f.close()

    def add_refuel_volume(self, iter):
        'Returns h, the total fuel salt height in the gas plenum after refuel salt has been added.'
        'The V0 variable is the initial volume from the BOC core.'
        filename = "Salt_volume_dep_step_{}.out".format(iter)
        volume = genfromtxt("{}".format(filename), delimiter='')   # Total fuel salt volume for current depletion step
        total_refuel_vol = volume - self.V0
        h = (total_refuel_vol)/(math.pi*200**2)
        return h
    
    ########################################
    #      Write Parallel KENO Decks
    ########################################

    def write_qsub_file(self, iter):
        'Writes a single qsub file outside of the refuel directories.'
        deck_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(deck_path)
        shell_content = '''
#!/bin/bash

#PBS -V
#PBS -q xeon
#PBS -l nodes=1:ppn=64
#PBS -l mem=200gb

cd $PBS_O_WORKDIR

hostname

module unload mpi
module load openmpi/2.1.6
module load scale/6.3.1-mpi
export DATA=/opt/scale6.3_data

export HDF5_USE_FILE_LOCKING=FALSE

scalerte -m -N 32 EIRENE.inp'''
        s = open(self.qsub_name, 'w')
        s.write(shell_content)
        s.close()

    def write_KENO_decks(self, iter):
        'Creates a directory for each refuel amount and writes a KENO deck there.'
        run_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))            # Make the depletion step directory
        os.mkdir(run_path)
#        os.chdir(run_path)
        for x in range(len(self.rvols)):
#            os.chdir('.')
            #########################################################################################
            #  Change directory to the location of the .f71 in order to read the .f71 file properly
            #########################################################################################
            if iter == 1:
                os.chdir(self.BOC_f71_path)
            else:
                last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter-1))
                os.chdir(last_dep_step_path)
            scale_fuel = self.write_scale_mat(self.rvols[x],iter)
            os.chdir(run_path)
            path_KENO = f'refuel_{self.rvols[x]:5.01f}'
            if not os.path.isdir(path_KENO):
                os.mkdir(path_KENO)
            os.chdir(path_KENO)
            self.write_fuel_salt_volume(self.rvols[x],iter, path_KENO)
            h = 440.5 + self.add_refuel_volume(iter)
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
 npg={self.npg} nsk={self.nsk} gen={self.gen} sig={self.sig}
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

            shell_content = '''
#!/bin/bash

#PBS -V
#PBS -q xeon
#PBS -l nodes=1:ppn=64

cd $PBS_O_WORKDIR

hostname

module unload mpi
module load openmpi/2.1.6
module load scale/6.3.1-mpi
export DATA=/opt/scale6.3_data

export HDF5_USE_FILE_LOCKING=FALSE

scalerte -m -N 32 EIRENE.inp'''
            q = open(self.qsub_name, 'w')
            q.write(shell_content)
            q.close()

            os.system('qsub runEIRENE-Scale.sh')  # Submit job

            os.chdir('..')

    def add_shell_permission(self, iter):
        'Gives permission to execute the shell script.'
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(main_path)
        os.system('chmod +x runEIRENE-Scale.sh')

    def write_conv_data(self, iter):
        'Writes the shell script convert-data.sh that extracts k-eff data from the finished KENO runs and prints them to file data-temp.out.'
        shell_main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        content = '''
#!/bin/bash

grep 'best estimate' refuel_*/EIRENE.out | sed -E -e s/^refuel.//g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g  > data-temp.out'''
        os.chdir(shell_main_path)
        f = open("convert-data.sh", 'w')
        f.write(content)
        f.close()

    def add_shell_permission_cd(self, iter):
      'Gives permissino to execute the convert-data.sh shell script.'
      main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
      os.chdir(main_path)
      os.system('chmod +x convert-data.sh')

    def convert_data(self, iter):
        'Submits the convert-data.sh shell script to extract KENO k-eff data and print to data-temp.out file.'
        os.system('/home/cmoss9/EIRENE11/SD7pct/dep_step_{}/./convert-data.sh'.format(iter))
    
    def run_SCALE_KENO(self, iter):
        'Runs SCALE from laptop terminal.'
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(main_path)
        for x in range(len(self.rvols)):
            deckpath = f'refuel_{self.rvols[x]:5.01f}'
            os.chdir(deckpath)
            os.system('/home/cmoss9/EIRENE11/SD7pct/dep_step_{}/./runEIRENE-Scale.sh'.format(iter))
            os.chdir('..')

    def read_outfile(self, iter):
        'Checks whether the run is finished and reads k-eff data.'
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        check = False       # Initialize check
        while check == False:
          if os.path.exists(main_path+'/data-temp.out'):
            check = True
            break
          else:                 
              check = False
              print("The atoms are still working; please stand by...")
              time.sleep(10.0)
        filename = "data-temp.out"
        data = genfromtxt("{}".format(filename), delimiter='')
        return data
    
    def get_crit_refuel(self, iter):
        'Fits the k-eff data and returns the critical refuel amount.'
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(main_path)
        data = self.read_outfile(iter)
        refuel = data[:,0]     # Refuel amounts in cm3
        keff = data[:,1]       # k-eff data
        kerr = data[:,2]       # k-eff error
        # Fit the data:
        m, b = np.polyfit(refuel, keff, 1)
        self.crit_refuel = (1 - b)/m     # Refuel amount for which k-eff = 1
        
        ##### Check if critical refuel amount is negative #####

        if self.crit_refuel < 1:
            self.crit_refuel = 0
        else:
            pass        

        f = open("crit_refuel_{}.out".format(iter), 'w')
        f.write("{}".format(self.crit_refuel))
        f.close()
        
        return self.crit_refuel
    
        ##################################################################################################
        #   Write New TRITON Deck With Crit Refuel (May want to put this in a different module script)
        ##################################################################################################
    
    def write_new_power_dens(self,iter):
        'Uses the normalization factor from the previous depletion step to convert the normalized MTHM mass to true MTiHM mass, then calculates the new power density.'
        path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(path)
        filename = "MTiHM_dep_step_{}.out".format(iter)
        MTiHM = genfromtxt("{}".format(filename), delimiter='')   # New MTiHM value calculated and written to .out file for current depletion step
        new_power_dens = self.power/MTiHM
        return new_power_dens
    
    def write_new_TRITON_deck(self, iter):
        'Writes a new TRITON deck for the current depletion step AND a KENO deck for the end-of-cycle depletion.'
        rvol = self.get_crit_refuel(iter)
        if iter == 1:
            os.chdir(self.BOC_f71_path)
        else:
            last_dep_step = iter - 1
            last_dep_step_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(last_dep_step))
            os.chdir(last_dep_step_path)
    
        new_scale_fuel = self.write_scale_mat(rvol, iter)
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(main_path)
        self.write_fuel_salt_volume(rvol,iter, main_path)
        h = 440.5 + self.add_refuel_volume(iter)
        # Save Height to .out File
        f = open("TRITON_height_step_{}.out".format(iter), 'w')
        f.write("{}".format(h))
        f.close()

        self.write_MTiHM_file(iter)
        
        new_power_dens = self.write_new_power_dens(iter)    # New power density in MW/MTiHM
        new_triton_deck = f'''=t6-depl parm=(addnux=4)
EIRENE SCALE/TRITON model, UF4 mol% = {self.UF4molpct}
ce_v7.1

read comp
{new_scale_fuel}

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

' Burnup values. Power in MW/MTiHM
read burndata
   power={new_power_dens} burn=7 nlib=1 end
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
 npg=20000 nsk=20 gen=5020 sig=50e-5
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
end model

end
'''
          
        main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        try:
#            os.makedirs(self.deck_path, exist_ok=True)
            os.chdir(main_path)
            fh = open(self.deck_name, 'w')
            fh.write(new_triton_deck)
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ")
            print(e)

    def write_run_TRITON(self, iter):
        'Writes the shell script runEIRENE-SCALE.sh to run new KENO deck.'
        shell_main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        shell_content = '''
#!/bin/bash

scalerte EIRENE.inp'''
        try:                # Write the deck
            os.chdir(shell_main_path)
            f = open("runEIRENE-SCALE.sh", 'w')
            f.write(shell_content)
            f.close()
        except IOError as e:
            print("Unable to write to file", f)
            print(e)

    def add_shell_permission_run_TRITON(self, iter):
      'Gives permission to execute the convert-data.sh shell script.'
      shell_main_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
      os.chdir(shell_main_path)
      os.system('chmod +x runEIRENE-SCALE.sh')

    def run_SCALE(self, iter):
      'Submits SCALE job.'
      run_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))            # Make the depletion step directory
      os.chdir(run_path)
      os.system('qsub runEIRENE-Scale.sh')

###################################################################################################################################

if __name__ == '__main__':
    print("This module generates and runs parallel KENO decks for different refuel volume additions.")
    input("Press Ctrl+C to quit, or press enter to test it.")
    iter = int(input("Enter the depletion step : "))
    decks = Refuel_Deck()
    print("***** Writing SCALE KENO decks...")
    decks.write_KENO_decks(iter)

###################################################################################################################################


