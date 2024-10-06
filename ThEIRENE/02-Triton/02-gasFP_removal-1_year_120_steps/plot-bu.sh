#!/bin/bash

cat ThEIRENE.out  | awk '/best est/{print $6" "$10}' > _dat_keff
cat ThEIRENE.out  | awk '/Transport k=/{print $3" "$6" "$10" "$14}' > _dat_bu
paste _dat_bu _dat_keff > data_bu.txt
rm _dat_keff _dat_bu

gnuplot << EOGA
set term postscript portrait enhanced color size 15cm,12cm
set out "keff-days.ps"
unset bars
#set autoscale extend
set grid
set key bottom
set format x "%.0f"
set format y "%.2f"
set title "ThEIRENE FLiBe-U 5mol\% UF_4"
set xlabel "Effective full power days at 400 MW"
#set ylabel "ThEIRENE {/Symbol r} [pcm]"
set ylabel "ThEIRENE k_{eff}"

# Make sure X axis extends beyond the data points
stats 'data_bu.txt' u 1 name 'Xax' nooutput
Xax_step = (Xax_max-Xax_min) / Xax_records
Xmin = Xax_min - Xax_step
Xmax = Xax_max + Xax_step

#f(x) = a*x +b
#fit [:] f(x) 'data_bu.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5)  yerrors via a,b

#plot [Xmin:Xmax] 'data_bu.txt' u 1:(\$2-1)*1e5/\$2:(\$3*1e5) w e ls 7 ps .5 lc "blue" notit, [:] f(x) ls 1 lw 0.4 lc "navy" tit sprintf(" slope = %5.2f {/Symbol \261} %5.2f pcm/K", a, a_err)
plot [Xmin:Xmax] 'data_bu.txt' u 1:5:6 w e ls 7 ps .5 lc "red" notit

set out
EOGA

convert           \
   -verbose       \
   -density 500   \
   -trim          \
    keff-days.ps      \
   -quality 100   \
   -flatten       \
   -sharpen 0x1.0 \
   -geometry 1600x1000 \
    keff-days.png


