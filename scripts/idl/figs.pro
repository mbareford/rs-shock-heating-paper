@Start.pro

data_path = '../../data'
fig_path = '../../images/idl'


get_en_plot_data,data_path+'/128x128x256',100
get_en_plot_data,data_path+'/256x256x512',100
get_en_plot_data,data_path+'/512x512x1024',100

get_htvisc_data,data_path+'/256x256x512/midplanes',0,400,10, 0,256,0,256,256,256

do_rad_line_plot, data_path+'/256x256x512/midplanes',95, 0, 9,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 2,0.0d0,0, 27.0d0,0.5d0, 1.0d0,'r', '', data_path+'/256x256x512/midplanes/bz_r_95.txt'
do_rad_line_plot, data_path+'/256x256x512/midplanes',95, 0, 16,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 2,0.0d0,0, 27.0d0,0.5d0, 1.0d0,'r', '', data_path+'/256x256x512/midplanes/btheta_r_95.txt'


do_2d_colour_plot, data_path+'/256x256x512/midplanes',95, 0, 12,16,0,0.0d0,0.04d0,0,0,0, -1.0d0,1.0d0, -1.0d0,1.0d0, 2,0.0d0,0, '','','','', 1000.0d0,[190.0d0,50.0d0,130.0d0,80.0d0], {YTicks:5, YTickname:['-1', '-0.5', '0', '0.5', '1']}, fig_path+'/256x256x512/xy_visc_sca.ps'

do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 12,16,0,0.0d0,0.0d0,0,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','','','', 1000.0d0,[190.0d0,50.0d0,130.0d0,80.0d0], fig_path+'/256x256x512/visc_sca_rf_vec.ps'

do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 32,22,1,-0.11d0,0.11d0,1,0, 2,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/vz_sca_v_vec.ps'

do_2d_colour_plot_with_contours, data_path+'/256x256x512/midplanes',95, 0, 54,13,0,0.0d0,0.0d0,1,0, 12,'black',2.0d0,12,0.0d0,0.0d0, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/slow_mach_sca_visc_con.ps'

do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 33,22,1,-15.685d0,15.685d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/jz_sca_rf_vec.ps'

do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 66,16,0,0.0d0,0.02d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/plasmab_sca_rf_vec.ps'



do_2d_colour_plot, data_path+'/256x256x512',95, 0, 3,16,0,0.0d0,25.0d0,1,0,1, -1.0d0,1.0d0, -5.0d0,5.0d0, 1,0.0d0,0, '','',' ','', 1000.0d0,[150.0d0,100.0d0,400.0d0,100.0d0], {YTicks:2, YTickname:['-1', '0', '1'], XTicks:6, XTickname:['-5', '-3', '-1', '0', '1', '3', '5']}, fig_path+'/256x256x512/xz_temp_sca_95.ps'

do_2d_colour_plot, data_path+'/256x256x512',120, 0, 3,16,0,0.0d0,25.0d0,0,0,1, -1.0d0,1.0d0, -5.0d0,5.0d0, 1,0.0d0,0, '','','','', 1000.0d0,[150.0d0,100.0d0,400.0d0,100.0d0], {YTicks:2, YTickname:['-1', '0', '1'], XTicks:6, XTickname:['-5', '-3', '-1', '0', '1', '3', '5']}, fig_path+'/256x256x512/xz_temp_sca_120.ps'

do_2d_colour_plot, data_path+'/256x256x512',400, 0, 3,16,0,0.0d0,25.0d0,0,0,1, -1.0d0,1.0d0, -5.0d0,5.0d0, 1,0.0d0,0, '','','','', 1000.0d0,[150.0d0,100.0d0,400.0d0,100.0d0], {YTicks:2, YTickname:['-1', '0', '1'], XTicks:6, XTickname:['-5', '-3', '-1', '0', '1', '3', '5']}, fig_path+'/256x256x512/xz_temp_sca_400.ps'



;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 1,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/rho_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 3,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/temp_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 21,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/bmag_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 9,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/jmag_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 8,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/jcrit_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 10,16,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/jsuper_sca_rf_vec.ps'

;do_2d_colour_plot_with_vectors, data_path+'/256x256x512/midplanes',95, 0, 35,22,0,0.0d0,0.0d0,1,0, 5,0,1l,0.05d0,1, 0.0d0,0.81d0, 0.0d0,0.81d0, 2,0.0d0, '','',' ','', 1000.0d0,[190.0d0,50.0d0,130.0d0,340.0d0], fig_path+'/256x256x512/rfz_sca_rf_vec.ps'
