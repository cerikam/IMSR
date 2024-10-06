#!/bin/sh
#
# Plots k_eff(enrichment)
#
# Ondrej Chvala <ochvala@utk.edu>
# MIT license


grep 'best estimate' Upct_*/ThEIRENE.out | sed -E -e s/^Upct.//g -e 's/.ThEIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g | tee data_crit_ufrac.txt


gnuplot << EOGA
set term postscript portrait enhanced color size 15cm,12cm
set out "rho-crit_ufrac.ps"
unset bars
set autoscale extend
set grid
set key bottom
set format x "%.2f"
set format y "%.0f"
set title "ThEIRENE FLiBe-U-Th, 8mol\% (U+Th)F4, 19.75\% LEU\nSCALE/KENO"
set xlabel "Uranium mol %"
set ylabel "ThEIRENE {/Symbol r} [pcm]"

# Make sure X axis extends beyond the data points
stats 'data_crit_ufrac.txt' u 1 name 'Xax' nooutput
Xax_step = (Xax_max-Xax_min) / Xax_records
Xmin = Xax_min - Xax_step
Xmax = Xax_max + Xax_step

f(x) = a*x +b
fit [:] f(x) 'data_crit_ufrac.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5)  yerrors via a,b

#plot [Xmin:Xmax] 'data_crit_ufrac.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit, f(x) ls 1 lw 0.4 lc "navy" tit sprintf(" slope = %5.2f {/Symbol \261} %5.2f pcm/K", a, a_err)
#plot [Xmin:Xmax] 'data_crit_ufrac.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit, [:] f(x) ls 1 lw 0.4 lc "navy" tit sprintf(" slope = %5.0f {/Symbol \261} %5.0f pcm/Uenr\%", a, a_err)
plot [Xmin:Xmax] 'data_crit_ufrac.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit

set out
EOGA

convert           \
   -verbose       \
   -density 500   \
   -trim          \
    rho-crit_ufrac.ps      \
   -quality 100   \
   -flatten       \
   -sharpen 0x1.0 \
   -geometry 1600x1000 \
    rho-crit_ufrac.png



