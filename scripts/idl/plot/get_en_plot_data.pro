PRO get_en_plot_data, input_path, skip_i

en = getenergy(wkdir=input_path)
t_i_max = SIZE(en.time,/N_ELEMENTS) - 1

en_tot0 = en.en_b[0] + en.en_ke[0] + en.en_int[0]

lun = 10
OPENW, lun, input_path+"/en.txt", WIDTH=512
FOR t_i=LONG(0),t_i_max,skip_i DO BEGIN  
  
  en_tot = en.en_b[t_i] + en.en_ke[t_i] + en.en_int[t_i]
  en_dis = en_tot/en_tot0 - 1.0d0 

  PRINTF, lun, en.time[t_i], ' ', en.en_b[t_i], ' ', en.en_ke[t_i], ' ', en.en_int[t_i], ' ', $
          en.heating_visc[t_i], ' ', en.heating_ohmic[t_i], ' ', en.heating_dp[t_i], ' ', $ 
          en.con_supp[t_i], ' ',  en.x_loss_con[t_i], ' ', en.y_loss_con[t_i], ' ', en.z_loss_con[t_i], ' ', en.loss_rad[t_i], ' ', $
          en.max_jmag[t_i], ' ', en.max_temp[t_i], ' ', en.eta_crit_frac[t_i], ' ', $  
          en.htvisc_xy[t_i], ' ', en.htvisc_xz[t_i], ' ', en.htvisc_yz[t_i], ' ', $
          en.htvisc_xx[t_i], ' ', en.htvisc_yy[t_i], ' ', en.htvisc_zz[t_i], ' ', $
          en.rke_neg[t_i], ' ', en.rke_pos[t_i], ' ', en_dis
          
ENDFOR
CLOSE, lun

END
