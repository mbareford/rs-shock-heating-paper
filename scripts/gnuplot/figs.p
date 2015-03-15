set xrange [*:*]
set yrange [*:*]
set palette rgbformulae 33,13,10
set pm3d map
unset pm3d
unset ztics
unset key
unset label

mu0 = 4.0e-7*pi
mp = 1.6726e-27
me = 9.1094e-31
ep = 8.8542e-12
ec = 1.6022e-19
kB = 1.3807e-23
gamma = 5.0/3.0  
mf = 1.0
m0 = mf*mp

B0 = 5.0e-3 
L0 = 1.0e6  
n0 = 1.0e15 
    
rho0 = m0*n0     
vA = B0/sqrt(mu0*rho0) 
tA = L0/vA             

P0 = B0**2/mu0    
En0 = P0/rho0      
Tp0 = En0*(m0/kB)  
j0 = B0/(mu0*L0)  


data_path = '"../../data'
fig_path = '"../../images/gnuplot'


set macros

set style line 1 lt 1 lw 1
set style line 2 lt 2 lw 1
set style line 3 lt 3 lw 1
set style line 4 lt 4 lw 1
set style line 5 lt 5 lw 1
set style line 6 lt 6 lw 1
set style line 7 lt 7 lw 1
set style line 8 lt 8 lw 1
set style line 9 lt 1 lw 1 lc 9

set terminal postscript eps color enhanced 28 dl 4


lm = 1.8

set xrange [0:1]
set yrange [0:1]
set xlabel "{/Italic r}" font "Helvetica,30"
set label "{/Italic B}_{/Italic z}" font "Helvetica,30" at 0.8,0.9
set label "{/Italic B}_{/Symbol-Oblique q}" font "Helvetica,30" at 0.8,0.2
set output @fig_path/b.eps"
plot lm*x*((1.0-x**2)**3) with lines lt -1 lw 1, sqrt(((lm**2)/7.0)*((1.0-x**2)**7 - 1.0) - (lm**2)*(x**2)*(1.0-x**2)**6 + 1.0) with lines lt -1 lw 1
unset label

set xrange [0:1]
set yrange [*:*]
set xlabel "{/Italic r}" font "Helvetica,30"
set label "{/Symbol-Oblique a}" font "Helvetica,30" at 0.8,3.5
set output @fig_path/alpha.eps"
plot (1.0/sqrt(((lm**2)/7.0)*((1.0-x**2)**7 - 1.0) - (lm**2)*(x**2)*(1.0-x**2)**6 + 1.0))*2.0*lm*((1.0-x**2)**2)*(1.0-4.0*(x**2)) with lines lt -1 lw 1, 0 lt 9
unset label


unset xlabel
set xrange [0:400]
set yrange [*:*]
set ylabel "{/Italic W}" font "Helvetica,30"
set output @fig_path/en_w.eps"
set key at 375.0,87.7
set key spacing 2.0
set key font ",20"
plot @data_path/128x128x256/en.txt" using 1:2 with lines ls 1 lc 3 lw 2 title "128x128x256", @data_path/256x256x512/en.txt" using 1:2 with lines ls 1 lc -1 title "256x256x512", @data_path/512x512x1024/en.txt" using 1:2 with lines ls 1 lc 1 lw 2 title "512x512x1024", "./tags.txt" using ($1==1&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==1&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label
unset key

set yrange [*:*]
set ylabel "{/Italic E}_{kin}" font "Helvetica,30"
set output @fig_path/en_kin.eps"
plot @data_path/128x128x256/en.txt" using 1:3 with lines ls 1 lc 3 lw 2, @data_path/256x256x512/en.txt" using 1:3 with lines ls 1 lc -1, @data_path/512x512x1024/en.txt" using 1:3 with lines ls 1 lc 1 lw 2, "./tags.txt" using ($1==2&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==2&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label

set xlabel "{/Italic t}" font "Helvetica,30"
set yrange [*:*]
set ylabel "{/Italic E}_{int}" font "Helvetica,30"
set output @fig_path/en_int.eps"
plot @data_path/128x128x256/en.txt" using 1:4 with lines ls 1 lc 3 lw 2, @data_path/256x256x512/en.txt" using 1:4 with lines ls 1 lc -1, @data_path/512x512x1024/en.txt" using 1:4 with lines ls 1 lc 1 lw 2, "./tags.txt" using ($1==3&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==3&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label

set yrange [*:*]
set ylabel "{/Italic E}_{dis}" font "Helvetica,30"
set output @fig_path/en_dis.eps"
plot @data_path/128x128x256/en.txt" using 1:24 with lines ls 1 lc 3 lw 2, @data_path/256x256x512/en.txt" using 1:24 with lines ls 1 lc -1, @data_path/512x512x1024/en.txt" using 1:24 with lines ls 1 lc 1 lw 2, "./tags.txt" using ($1==4&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==4&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label

set yrange [*:*]
set ylabel "{/Italic Shock heating / Ohmic heating}" font "Helvetica,30"
set output @fig_path/ht_visc_ohmic.eps"
plot @data_path/128x128x256/en.txt" using 1:($6>0.0?$5/$6:1/0) with lines ls 1 lc 3 lw 2, @data_path/256x256x512/en.txt" using 1:($6>0.0?$5/$6:1/0) with lines ls 1 lc -1, @data_path/512x512x1024/en.txt" using 1:($6>0.0?$5/$6:1/0) with lines ls 1 lc 1 lw 2, "./tags.txt" using ($1==5&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==5&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label

set xlabel "{/Italic t}" font "Helvetica,30"
set xrange [*:*]
set yrange [*:*]
set ylabel "{/Italic j}_{max}" font "Helvetica,30"
set output @fig_path/j_max.eps"
plot @data_path/128x128x256/en.txt" using 1:13 with lines ls 1 lc 3 lw 2, @data_path/256x256x512/en.txt" using 1:13 with lines ls 1 lc -1, @data_path/512x512x1024/en.txt" using 1:13 with lines ls 1 lc 1 lw 2, "./tags.txt" using ($1==6&&$2==1?$3:1/0):4 with points ls 6 lc 9 lw 2 ps 2 title "", "./tags.txt" using ($1==6&&$2==3?$3:1/0):4 with points ls 8 lc 9 lw 2 ps 2 title ""
unset label

set xrange [0:400]
set yrange [*:*]
set key at 120,0.475
set key spacing 1.5
set xlabel "{/Italic t}" font "Helvetica,30"
set ylabel "{/Italic Shock Heating (tensor comps.)}" font "Helvetica,30"
set output @fig_path/256x256x512/ht_visc_comp.eps"
plot @data_path/256x256x512/en.txt" using 1:16 with lines ls 1 lc -1 lw 2 title "{/Italic xy}", @data_path/256x256x512/en.txt" using 1:19 with lines ls 6 lc -1 lw 2 title "{/Italic xx}", @data_path/256x256x512/en.txt" using 1:20 with lines ls 2 lc -1 lw 2 title "{/Italic yy}", @data_path/256x256x512/en.txt" using 1:17 with lines ls 3 lc -1 lw 2 title "{/Italic xz}", @data_path/256x256x512/en.txt" using 1:18 with lines ls 5 lc -1 lw 2 title "{/Italic yz}", @data_path/256x256x512/en.txt" using 1:21 with lines ls 4 lc -1 lw 2 title "{/Italic zz}"
unset label

unset ylabel
unset key
set label "{/Italic z} = 0" font "Helvetica,30" at 330.0,0.008
set output @fig_path/256x256x512/ht_visc_comp_apex.eps"
plot @data_path/256x256x512/midplanes/visc_apex.txt" using 1:5 with lines ls 1 lc -1 lw 2, @data_path/256x256x512/midplanes/visc_apex.txt" using 1:3 with lines ls 6 lc -1 lw 2, @data_path/256x256x512/midplanes/visc_apex.txt" using 1:4 with lines ls 2 lc -1 lw 2, @data_path/256x256x512/midplanes/visc_apex.txt" using 1:7 with lines ls 3 lc -1 lw 2, @data_path/256x256x512/midplanes/visc_apex.txt" using 1:6 with lines ls 5 lc -1 lw 2, @data_path/256x256x512/midplanes/visc_apex.txt" using 1:8 with lines ls 4 lc -1 lw 2
unset label
unset xlabel


set xrange [0.0:1.0]
set yrange [0.0:0.8]
set xlabel "{/Italic r}" font "Helvetica,30"
set ylabel "{/Italic B}_{/Symbol-Oblique q}" font "Helvetica,30"
set key at 0.35,0.75
set key spacing 2.0
set key font ",25"
set output @fig_path/256x256x512/btheta_r_95.eps"
plot lm*x*((1.0-x**2)**3) with lines ls 4 lc 9 lw 2 title "0 {/Italic t}_{A}", @data_path/256x256x512/midplanes/btheta_r_95.txt" using 1:2 with lines ls 1 lc 1 lw 2 title "95 {/Italic t}_{A}"
unset label
unset xlabel
unset title
unset key

set xrange [0.0:1.0]
set yrange [0.4:1.0]
set xlabel "{/Italic r}" font "Helvetica,30"
set ylabel "{/Italic B}_{/Italic z}" font "Helvetica,30"
set output @fig_path/256x256x512/bz_r_95.eps"
plot sqrt(((lm**2)/7.0)*((1.0-x**2)**7 - 1.0) - (lm**2)*(x**2)*(1.0-x**2)**6 + 1.0) with lines ls 4 lc 9 lw 2, @data_path/256x256x512/midplanes/bz_r_95.txt" using 1:2 with lines ls 1 lc 1 lw 2
unset label
unset xlabel
unset title
unset key
