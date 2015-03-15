; get the volume-integrated viscous heating over a subsection of the grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO get_htvisc_data, energy_input_path, input_path, ss_i_min, ss_i_max, ss_i_step, x_i_min, x_i_max, y_i_min, y_i_max, z_i_min, z_i_max, dV
        
  en = getenergy(wkdir=energy_input_path)
  t_i_max = SIZE(en.time,/N_ELEMENTS) - 1

  lun = 10
  OPENW, lun, input_path+'/visc_apex.txt', WIDTH=256  
  FOR ss_i = ss_i_min,ss_i_max,ss_i_step DO BEGIN
  
    PRINT, 'ss_i=', ss_i 

    t_i = 1
    REPEAT BEGIN
      dt = en.time[t_i] - en.time[t_i-1]
      t_i = t_i + 1
    ENDREP UNTIL en.time[t_i] GE ss_i || t_i GT t_i_max
    

    ds = getdata(ss_i, wkdir=input_path)
    
    visc_heat_xx = TOTAL(ds.visc_heat_xx[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    visc_heat_yy = TOTAL(ds.visc_heat_yy[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    visc_heat_xy = TOTAL(ds.visc_heat_xy[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    visc_heat_xz = TOTAL(ds.visc_heat_xz[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    visc_heat_yz = TOTAL(ds.visc_heat_yz[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    visc_heat_zz = TOTAL(ds.visc_heat_zz[x_i_min:x_i_max,y_i_min:y_i_max,z_i_min:z_i_max])*dV*dt
    
    visc_heat = visc_heat_xx + visc_heat_yy + visc_heat_xy + visc_heat_xz + visc_heat_yz + visc_heat_zz
    
    PRINTF, lun, ss_i, ' ', visc_heat, ' ', visc_heat_xx, ' ', visc_heat_yy, ' ', visc_heat_xy, ' ', visc_heat_xz, ' ', visc_heat_yz, ' ', visc_heat_zz
      
  ENDFOR
  CLOSE, lun
        
END

