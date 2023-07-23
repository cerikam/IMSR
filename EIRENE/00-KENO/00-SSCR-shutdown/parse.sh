#!/bin/bash

grep 'best estimate' -- */*//EIRENE.out | sed -E -e s/.enr./\ \ /g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g > data.txt


# 01-SS-up 2.650 1.00087 0.00049
# 02-SS-mid 2.650 0.97811 0.00049

# Difference between SSCRs up/middle is (86.92437579305269 - 2237.9895921726556) = 2325 +- 70 pcm
