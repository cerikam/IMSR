#!/bin/env python

"""
Utilities for EIRENE in SCALE project
Ondrej Chvala <ochvala@utexas.edu>
"""

import re
debug:int = 0

def parse_mixing_table(filename):
    """Parses input file and extracts mixing table, and total volume and mass table"""
    table_start_pattern = re.compile(r'mixture\s*=\s*(\d+)')
    column_names = []
    mixing_table = {}
    current_mixture = None

    file = open(filename, 'r')
    # Parse Mixing table
    for line in file:
        line = line.strip()
        if line == "":
            continue

        if current_mixture is not None:
            if line.find('finished preparing the keno input data') > 0 or line.find('=================') > -1:
                if debug > 5:
                    print("DONE:", line)
                break

        # Check if the line starts a new mixture
        match = table_start_pattern.match(line)
        if match:
            if debug > 5:
                print(line)
            mixture_number = int(match.group(1))
            current_mixture = mixing_table[mixture_number] = {}
            continue

        # Check if the line contains column names
        if line.startswith('nuclide'):
            line = line.replace('wgt. frac.', 'wgt-frac.')
            line = line.replace('nuclide', '')
            column_names = line.split()
            continue

        # Parse the data rows
        if current_mixture is not None and column_names:
            values = line.split()
            nuclide = values.pop(0)
            current_mixture[nuclide] = {column: value for column, value in zip(column_names, values)}

    # Parse volumes and masses
    current_mixture = None
    file.seek(0)
    for line in file:
        line = line.strip()
        if line == "":
            continue

        if current_mixture is not None:
            if line.find('-------') > -1:
                if debug > 2:
                    print("DONE:", line)
                break
            values = line.split()
            mixture = int(values[0])
            mixing_table[mixture]['tot_volume'] = values[1]
            mixing_table[mixture]['tot_volume_err'] = values[3]
            mixing_table[mixture]['tot_mass'] = values[4]
            mixing_table[mixture]['tot_mass_err'] = values[6]

        if line.find('total mixture volume ') > 0:
            current_mixture = "table found"

    return mixing_table


def get_MTiHM(mixing_table):
    """Returns Metric Tonnes of Initial Heavy Metal for actinide-bearning mixtures"""
    MTiHM = {}
    for m in mixing_table.keys():
        nuclides = [x for x in mixing_table[m].keys() if x.isnumeric()]
        HMs = [x for x in nuclides if 1000000 > int(x) >= 90000]
        if HMs:
            MTiHM[m] = 0.0
            for HM in HMs:
                tot_mix_mass = float(mixing_table[m]['tot_mass'])
                nuc_wgt_frac = float(mixing_table[m][HM]['wgt-frac.'])
                MTiHM[m] += 1e-6 * tot_mix_mass * nuc_wgt_frac

    return MTiHM


if __name__ == '__main__':
    filename = 'IMSRRev6.out'
    mixing_table = parse_mixing_table(filename)
    print(mixing_table)
    mtihm = get_MTiHM(mixing_table)
    print(mtihm)
