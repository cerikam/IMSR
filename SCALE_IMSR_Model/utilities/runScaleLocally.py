#!/bin/env python3

"""Runs all cases in the current directory locally
   Note, this runs serial SCALE! """

import os

SCALE_bin_path = os.getenv('SCALE_BIN', '/opt/scale/scale6.3.1/bin/')
if SCALE_bin_path == '/opt/scale/scale6.3.1/bin/':
    print(f"""WARNING: using default SCALE bin path:
         SCALE_BIN = {SCALE_bin_path}
""")


def main():
    # Stating directory
    cwd = os.getcwd()

    # Loop through all directories and sub-directories recursively
    for dir1, subdirs, files in os.walk(cwd):
        # Check if "SCALE_FILE" subdirectory exists in the current directory
        if "SCALE_FILE" in subdirs:
            os.chdir(f'{dir1}/SCALE_FILE/')
            print(f'Running {dir1}')
            os.system(f'{SCALE_bin_path}/scalerte -m SCALE_FILE.inp &')
            os.chdir(f'{cwd}')


if __name__ == '__main__':
    main()
