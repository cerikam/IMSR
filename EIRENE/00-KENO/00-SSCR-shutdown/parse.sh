#!/bin/bash

grep 'best estimate' -- */*//EIRENE.out | sed -E -e s/.enr./\ \ /g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g > data.txt


# 01-SS-up 2.650 1.00180 0.00049
# 02-SS-mid 2.650 0.97859 0.00049

# Difference between SSCRs up/middle is (179.67658215212856 + 2187.8416905956574) = 2367.518272747786 pcm
