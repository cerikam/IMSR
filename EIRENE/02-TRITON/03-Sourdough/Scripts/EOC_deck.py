# ***************************************************************************************************
#
#                        EIRENE Sourdough Write and Run KENO for EOC Depletion
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
import math
import sys, re
import subprocess

SCALE_bin_path: str = os.getenv('SCALE_BIN', '/opt/scale6.3.1/bin/')


#####################################################
#           Initialize EOC KENO Deck
#####################################################

# Initialize and create the input deck

class EOC_Deck(object):
    'Initial EIRENE EOC KENO deck.'
    def __init__(self):
        ##################################
        # EOC core and running parameters
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
        self.renrich = 0.07                                    # Refuel enrichment fraction

        # ************************
        #    Set KENO Params.
        # ************************
        self.npg:float = 20000           # Number per generation (npg)
        self.nsk:float = 20              # Number of skipped generations (nsk)
        self.gen:float = 5020            # Number of generations (gen)
        self.sig:float = 150e-5           # (sig)
        self.f71_name:str = 'EIRENE.f71'            # Name of .f71 TRITON file to read
        self.init_MTiHM:float = 6.883864319048779        # Initial MTiHM for BOC core
        
        # ************************
        #    Running Parameters
        # ************************
        self.queue:str = 'fill'                     # NECluster queue
        self.ppn:int = 64                           # ppn core count
        self.deck_name:str = 'EIRENE_EOC.inp'           # TRITON input file name
        self.deck_path:str = os.path.expanduser('~/EIRENE11/SD7pct/BOC')                     # Where to run the TRITON deck
        self.qsub_path:str = os.path.expanduser('~/EIRENE11/SD7pct/BOC/runEIRENE-Scale.sh')  # Full path to the qsub script
        self.power:float = 400.0                    # Power in MW
        self.MTiHM:float = 6.883864319048779        # Initial MTiHM for BOC core
        self.qsub_name:str = 'runEIRENE-Scale.sh'


    def get_EOC_burned_salt_adens(self):
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
        
        
        return burned_salt_vector

    def read_TRITON_height(self, iter):
        'Reads the refuel height file for the current depletion step (for depletion steps past the BOC run).'
        filename = "TRITON_height_step_{}.out".format(iter)
        salt_plenum_height = genfromtxt("{}".format(filename), delimiter='')
        return salt_plenum_height

    def write_EOC_scale_mat(self):
        'Returns a SCALE material composition block for the refuel salt mixed in with the burned salt.'
        enrich_percent = self.renrich*100    # Enrichment percent for refuel
        burned_salt_adens = self.get_EOC_burned_salt_adens()
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
            scale_mat += "{:10s} 1 0  {:>5e} 923.15 end \n".format(isotopes[i], burned_salt_adens[i])
        return scale_mat
    
    def write_EOC_KENO(self, iter):
        'Writes a new KENO deck for the EOC depletion.'
        run_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(run_path)
        scale_fuel = self.write_EOC_scale_mat()
        h = self.read_TRITON_height(iter)
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

    def run_SCALE(self, iter):
      'Submits SCALE job.'
      run_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))            # Make the depletion step directory
      os.chdir(run_path)
      os.system('qsub runEIRENE-Scale.sh')


    def write_qsub_file(self, iter):
        'Writes a single qsub file outside of the refuel directories.'
        deck_path = os.path.expanduser('~/EIRENE11/SD7pct/dep_step_{}'.format(iter))
        os.chdir(deck_path)
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
        s = open(self.qsub_name, 'w')
        s.write(shell_content)
        s.close()

###############################################################################################

if __name__ == '__main__':
    print("This module generates and runs an EOC KENO deck for EIRENE sourdough depletion.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    iter = int(input("Enter the depletion step : "))
    eoc = EOC_Deck()
    print("***** Writing EOC KENO deck...")
    eoc.write_EOC_KENO(iter)
    print("Success!")

################################################################################################
