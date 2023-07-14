#!/bin/env python3


""" Serpent model of example Integral REactor for Nuclear Education (SIRENE)
Script to find critical enrichment of SIRENE using EIRENE's atom densities
Ondrej Chvala <ochvala@utk.edu>
MIT license
"""

import salts
import os
import numpy as np


def tempK(tempC: float) -> float:
    return tempC + 273.15


def tempC(tempK: float) -> float:
    return tempK - 273.15


