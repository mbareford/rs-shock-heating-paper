1: [left] ./images/gnuplot/b, [right] ./images/gnuplot/alpha

2: [top left] ./images/gnuplot/en_w, [top right] ./images/gnuplot/en_kin,
   [bottom left] ./images/gnuplot/en_int, [bottom right] ./images/gnuplot/en_dis}

3: [left] ./images/gnuplot/ht_visc_ohmic, [right] ./images/gnuplot/j_max

4: [left] ./images/gnuplot/256x256x512/ht_visc_comp, [right] ./images/gnuplot/256x256x512/ht_visc_comp_apex (*** need complete set of snapshot data, 0-400 tA ***)

5: [left] ./images/inkscape/256x256x512/xy_visc_sca, [right] ./images/inkscape/256x256x512/visc_sca_rf_vec

6: [left] ./images/inkscape/256x256x512/vz_sca_v_vec, [right] ./images/inkscape/256x256x512/slow_mach_sca_visc_con

7: [left] ./images/inkscape/256x256x512/plasmab_sca_rf_vec, [right] ./images/inkscape/256x256x512/jz_sca_rf_vec

8: [left] ./images/gnuplot/256x256x512/btheta_r_19, [right] ./images/gnuplot/256x256x512/bz_r_19

9: ./images/inkscape/256x256x512/xz_temp_sca

10: ./images/inkscape/256x256x512/bfield_voheat


Instructions
------------
cd ./code/lare3d
# the readme.txt contains brief instructions on how to compile the source code, and then to execute the lare3d
# simulation in order to generate the CFD files, which are too big to include in the github repository at
# https://github.com/mbareford/rs-shock-heating-paper

cd ./data/256x256x512
mkdir midplanes
# change the names of the CFD files, so that each snapshot name is the Alfven time at which the snapshot was taken
mv 0080.cfd 0400.cfd
mv 0079.cfd 0395.cfd
...
mv 0019.cfd 0095.cfd
...
mv 0002.cfd 0010.cfd
mv 0001.cfd 0005.cfd

cd ../../scripts/python
# extract a z slices (-1 <= z <= 1) from the snapshots taken by the 256x256x512 resolution simulation
python extract_midplanes.py  

cd ../idl/
# generate the txt files that will be used by a gnuplot script to produce the eps files for figures 2-4 and 8
# generate the ps files for figures 5-7 and 9
idl figs.pro

cd ../visit/
# use visit v1.12.2 to generate the JPEGs for figure 10, simply load the session file (visit.session)
# and navigate to the folder containing the non-sliced CFD files.
# note, the lare3d visit plugin must be installed before visit can read CFD files,
# see http://ccpforge.cse.rl.ac.uk/gf/project/lare3d for plugin source code.

cd ../../images/inkscape/256x256x512
# use inkscape v0.48.4 r9939 to generate the PDFs for figures 5-7 and 9-10
# inkscape files are *.svg files, the content is added by importing postscript files, see ./images/inkscape/idl/256x256x512 (figs 5-7,9),
# and by importing jpeg files, see ./images/inkscape/visit/256x256x512 (fig 10).

cd ../../scripts/gnuplot
gnuplot figs.p
