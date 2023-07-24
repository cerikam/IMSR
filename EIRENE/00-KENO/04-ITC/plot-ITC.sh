#!/bin/sh
#
# Plots k_eff(enrichment)
#
# Ondrej Chvala <ochvala@utk.edu>
# MIT license


grep 'best estimate' ITC_*/EIRENE.out | sed -E -e s/^ITC.//g -e 's/.EIRENE.*eff//g'  -e 's/\+ or \-//g' -e 's/\*\*\*//g' -e 's/\s+/ /g' | sort -g  > data-temp.out

gnuplot << EOGA
set term postscript portrait enhanced color size 15cm,12cm
set out "rho-temp-ITC.ps"
unset bars
set grid
set key bottom
set format x "%.0f"
set format y "%.0f"
set title "EIRENE FLiBe-U 5mol\% UF_4, SCALE/KENO"
set xlabel "Core temperature [°C]"
set ylabel "EIRENE {/Symbol r} [pcm]"

# Make sure X axis extends beyond the data points
stats 'data-temp.out' u 1 name 'Xax' nooutput
Xax_step = (Xax_max-Xax_min) / Xax_records
Xmin = Xax_min - Xax_step
Xmax = Xax_max + Xax_step

f(x) = a*x +b
fit [:] f(x) 'data-temp.out' u 1:(\$2-1)*1e5/\$2:(\$3*1e5)  yerrors via a,b

#plot [Xmin:Xmax] 'data-temp.out' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit, f(x) ls 1 lw 0.4 lc "navy" tit sprintf(" slope = %5.2f {/Symbol \261} %5.2f pcm/K", a, a_err)
plot [Xmin:Xmax] 'data-temp.out' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit, [:] f(x) ls 1 lw 0.4 lc "navy" tit sprintf(" slope = %5.2f {/Symbol \261} %5.2f pcm/K", a, a_err)

set out
EOGA

convert           \
   -verbose       \
   -density 500   \
   -trim          \
    rho-temp-ITC.ps      \
   -quality 100   \
   -flatten       \
   -sharpen 0x1.0 \
   -geometry 1600x1000 \
    rho-temp-ITC.png



