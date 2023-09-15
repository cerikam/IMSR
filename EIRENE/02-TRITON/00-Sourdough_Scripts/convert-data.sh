#!/bin/sh
#
#

grep 'best estimate' refuel_*/EIRENE.out | sed -E -e s/^refuel.//g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g  > data-temp.out
