#!/bin/env python3

"""Runs all cases in the current directory on a PBS cluster"""

import os


def main():
    # Stating directory
    cwd = os.getcwd()

    # Loop through all directories and sub-directories recursively
    for dir1, subdirs, files in os.walk(cwd):
        # Check if "SCALE_FILE" subdirectory exists in the current directory
        if "SCALE_FILE" in subdirs:
            os.chdir(f'{dir1}/SCALE_FILE/')
            print(f'Running {dir1}')
            os.system(f'qsub {cwd}/runScale.sh')
            os.chdir(f'{cwd}')


if __name__ == '__main__':
    main()
