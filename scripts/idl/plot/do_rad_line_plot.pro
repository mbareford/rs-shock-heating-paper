FUNCTION index, ds_ext, d_i, c    
  c_arr = 0l
           
  CASE d_i OF
    0: BEGIN
         c_arr = ds_ext.x
       END
    1: BEGIN         
         c_arr = ds_ext.y
       END
    ELSE: BEGIN
            c_arr = ds_ext.z
          END
  ENDCASE
    
  c_i_min = 0l
  c_i_max = SIZE(c_arr, /N_ELEMENTS) - 1l
  
  c_i = c_i_min
  IF (c GE c_arr[c_i_min]) THEN BEGIN
    c_i = c_i_max
    IF (c LE c_arr[c_i_max]) THEN BEGIN
    
      cnt = 0l
      ge_c = WHERE(c_arr GE c, cnt)
      IF (cnt GT 0l) THEN BEGIN
        c_i = ge_c[0l]      
        IF (c_i GT c_i_min) THEN BEGIN
          IF (ABS(c-c_arr[c_i-1l]) LT ABS(c-c_arr[c_i])) THEN BEGIN  
            c_i = c_i - 1l
          ENDIF
        ENDIF
      ENDIF
      
    ENDIF
  ENDIF

  RETURN, c_i
END




;do_rad_line_plot, data_path+'/job03/256x256x512',19, 0, 9,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 2,0.0d0,0, 0.0d0,1.0d0, 1.0d0,'r', '', data_path+'/job03/256x256x512/bz_r_19.txt'
PRO do_rad_line_plot, input_path,t_i, cen, sca_id,sca_min,sca_max, d1_min,d1_max,d1_org, d2_min,d2_max,d2_org, d3_id,d3,d3_sum, theta_sel,theta_tol, rad_max,rad_title, output_filename, data_filename

  ds = getdata(t_i, wkdir=input_path,/rho,/bz)
  ;PLOT, ds.grid.x[128:256], ds.bz[128:256,128,256], psym=-2, xr=[0,1]
  
  IF (cen EQ 0) THEN BEGIN
    ds_ext = get_ext_ccdata(ds,/rho)
  ENDIF ELSE BEGIN
    ds_ext = get_ext_vcdata(ds,/rho)
  ENDELSE  
  
  d3_i = index(ds_ext,d3_id,d3)
  d3_i_min = d3_i
  d3_i_max = d3_i
  IF d3_sum EQ 1 THEN BEGIN
    d3_i_min = 0l
  ENDIF
  
  
  CASE d3_id OF
    0: BEGIN         
         d1_i_min = index(ds_ext,1,d1_min)
         d1_i_max = index(ds_ext,1,d1_max)
         d1_i_org = index(ds_ext,1,d1_org)
         d2_i_min = index(ds_ext,2,d2_min)
         d2_i_max = index(ds_ext,2,d2_max)
         d2_i_org = index(ds_ext,2,d2_org)
         IF d3_sum EQ 1 THEN BEGIN
           d3_i_max = SIZE(ds_ext.x, /N_ELEMENTS) - 1l
         ENDIF         
       END
    1: BEGIN         
         d1_i_min = index(ds_ext,0,d1_min)
         d1_i_max = index(ds_ext,0,d1_max)
         d1_i_org = index(ds_ext,0,d1_org)
         d2_i_min = index(ds_ext,2,d2_min)
         d2_i_max = index(ds_ext,2,d2_max)
         d2_i_org = index(ds_ext,2,d2_org)
         IF d3_sum EQ 1 THEN BEGIN
           d3_i_max = SIZE(ds_ext.y, /N_ELEMENTS) - 1l
         ENDIF
       END
    ELSE: BEGIN
            d1_i_min = index(ds_ext,0,d1_min)
            d1_i_max = index(ds_ext,0,d1_max)
            d1_i_org = index(ds_ext,0,d1_org)
            d2_i_min = index(ds_ext,1,d2_min)
            d2_i_max = index(ds_ext,1,d2_max)
            d2_i_org = index(ds_ext,1,d2_org)
            IF d3_sum EQ 1 THEN BEGIN
              d3_i_max = SIZE(ds_ext.z, /N_ELEMENTS) - 1l
            ENDIF
          END
  ENDCASE
     

       
  PRINT, 'd1_i_min=',d1_i_min
  PRINT, 'd1_i_max=',d1_i_max
  PRINT, 'd1_i_org=',d1_i_org
  PRINT, 'd2_i_min=',d2_i_min
  PRINT, 'd2_i_max=',d2_i_max
  PRINT, 'd2_i_org=',d2_i_org
  PRINT, 'd3_i_min=',d3_i_min
  PRINT, 'd3_i_max=',d3_i_max
  PRINT, 'd3_i=',d3_i
    
    
  ds = 0
  ds_ext = 0
  
  sca_3dv = 0    
  
  
  sca_2dp = DBLARR(d1_i_max-d1_i_min+1l,d2_i_max-d2_i_min+1l)
  theta_2dp = DBLARR(d1_i_max-d1_i_min+1l,d2_i_max-d2_i_min+1l)
  rad_2dp = DBLARR(d1_i_max-d1_i_min+1l,d2_i_max-d2_i_min+1l)
  
  
  CASE sca_id OF    
    1: BEGIN
         sca_title = 'density' 
         ds = getdata(t_i,wkdir=input_path,/rho)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho)
           sca_3dv = ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho)
           sca_3dv = ds_ext.rho
         ENDELSE                           
       END
    2: BEGIN
         sca_title = 'viscous heating'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat)        
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/visc_heat)
           sca_3dv = ds.visc_heat/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/visc_heat)
           sca_3dv = ds_ext.visc_heat/ds_ext.rho
         ENDELSE
       END
    3: BEGIN
         sca_title = 'viscous heating [xy]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xy_visc_heat)
           sca_3dv = ds.visc_heat_xy/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xy_visc_heat)
           sca_3dv = ds_ext.visc_heat_xy/ds_ext.rho
         ENDELSE
       END
    4: BEGIN
         sca_title = 'viscous heating [xz]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xz_visc_heat)
           sca_3dv = ds.visc_heat_xz/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xz_visc_heat)
           sca_3dv = ds_ext.visc_heat_xz/ds_ext.rho
         ENDELSE
       END
    5: BEGIN
         sca_title = 'viscous heating [yz]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_yz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/yz_visc_heat)
           sca_3dv = ds.visc_heat_yz/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/yz_visc_heat)
           sca_3dv = ds_ext.visc_heat_yz/ds_ext.rho
         ENDELSE
       END
    6: BEGIN
         sca_title = 'viscous heating [xx]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xx)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xx_visc_heat)
           sca_3dv = ds.visc_heat_xx/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xx_visc_heat)
           sca_3dv = ds_ext.visc_heat_xx/ds_ext.rho
         ENDELSE
       END
    7: BEGIN
         sca_title = 'viscous heating [yy]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_yy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/yy_visc_heat)
           sca_3dv = ds.visc_heat_yy/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/yy_visc_heat)
           sca_3dv = ds_ext.visc_heat_yy/ds_ext.rho
         ENDELSE
       END
    8: BEGIN
         sca_title = 'viscous heating [zz]'
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_zz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/zz_visc_heat)
           sca_3dv = ds.visc_heat_zz/ds.rho
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/zz_visc_heat)
           sca_3dv = ds_ext.visc_heat_zz/ds_ext.rho
         ENDELSE
       END
    9: BEGIN
         sca_title = 'bz'
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
           sca_3dv = ds_ext.bz
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)
           sca_3dv = ds_ext.bz
         ENDELSE
       END
   10: BEGIN
         sca_title = 'vz'
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity)
           sca_3dv = ds_ext.vz
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity)
           sca_3dv = ds.vz
         ENDELSE
       END
   11: BEGIN
         sca_title = 'jz'
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current)
           sca_3dv = ds_ext.jz
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current)
           sca_3dv = ds_ext.jz
         ENDELSE
       END
   12: BEGIN
         sca_title = 'fz'
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current,/lorentz)
           sca_3dv = ds_ext.fz
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current,/lorentz)
           sca_3dv = ds_ext.fz
         ENDELSE
       END
   13: BEGIN
         sca_title = 'cs mach number'
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/pressure,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/pressure,/velocity)
         ENDELSE
         sca_3dv = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)/SQRT(5.0d0/3.0d0*ds_ext.p/ds_ext.rho)
         nonshock_cnt = 0l
         nonshock_i = WHERE(sca_3dv LT 1.0d0, nonshock_cnt)
         IF (nonshock_cnt GT 0l) THEN BEGIN
           sca_3dv[nonshock_i] = 0.0d0
         ENDIF
       END
   14: BEGIN
         sca_title = 'v . f'
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity,/bfield,/current,/lorentz)
           sca_3dv = ds_ext.vx*ds_ext.fx + ds_ext.vy*ds_ext.fy + ds_ext.vz*ds_ext.fz
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity,/bfield,/current,/lorentz)
           sca_3dv = ds_ext.vx*ds_ext.fx + ds_ext.vy*ds_ext.fy + ds_ext.vz*ds_ext.fz
         ENDELSE
       END
   15: BEGIN
         sca_title = 'temperature'
         ds = getdata(t_i,wkdir=input_path,/temperature)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/temperature)
           sca_3dv = ds.temperature           
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/temperature)
           sca_3dv = ds_ext.temperature
         ENDELSE  
       END  
   16: BEGIN
         sca_title = 'btheta'
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)           
         ENDELSE
         
         nx = SIZE(ds.grid.x,/N_ELEMENTS) 
         ny = SIZE(ds.grid.y,/N_ELEMENTS)
         sca_3dv = ds_ext.bz
  
         FOR x_i = 0,nx-1l,1l DO BEGIN
           FOR y_i = 0,ny-1l,1l DO BEGIN
             x = ds_ext.x[x_i]
             y = ds_ext.y[y_i]      
         
             theta = (x EQ 0 AND y EQ 0) ? 0.0d0 : ATAN2(y,x)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360 - ABS(theta)
             ENDIF
             theta = theta*!dpi/180.0d0
         
             bx = ds_ext.bx[x_i,y_i,*]
             by = ds_ext.by[x_i,y_i,*]
             sca_3dv[x_i,y_i,*] = ds_ext.by[x_i,y_i,*]*COS(theta) - ds_ext.bx[x_i,y_i,*]*SIN(theta)
           ENDFOR        
         ENDFOR         
       END 
   17: BEGIN
         sca_title = 'br'
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)           
         ENDELSE
         
         nx = SIZE(ds.grid.x,/N_ELEMENTS) 
         ny = SIZE(ds.grid.y,/N_ELEMENTS)
         sca_3dv = ds_ext.bz
  
         FOR x_i = 0,nx-1l,1l DO BEGIN
           FOR y_i = 0,ny-1l,1l DO BEGIN
             x = ds_ext.x[x_i]
             y = ds_ext.y[y_i]      
         
             theta = (x EQ 0 AND y EQ 0) ? 0.0d0 : ATAN2(y,x)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360 - ABS(theta)
             ENDIF
             theta = theta*!dpi/180.0d0
         
             bx = ds_ext.bx[x_i,y_i,*]
             by = ds_ext.by[x_i,y_i,*]
             sca_3dv[x_i,y_i,*] = ds_ext.by[x_i,y_i,*]*COS(theta) + ds_ext.bx[x_i,y_i,*]*SIN(theta)
           ENDFOR        
         ENDFOR         
       END 
   ELSE: BEGIN
           PRINT, 'Error: unrecognised scalar id.'
           RETURN           
         END
  ENDCASE
       
  
  CASE d3_id OF
    0: BEGIN         
         d1_min = ds_ext.y[d1_i_min]
         d1_max = ds_ext.y[d1_i_max]
         d2_min = ds_ext.z[d2_i_min]
         d2_max = ds_ext.z[d2_i_max]
         dd1 = ABS(ds_ext.y[d1_i_min]-ds_ext.y[d1_i_min+1l])
         dd2 = ABS(ds_ext.z[d2_i_min]-ds_ext.z[d2_i_min+1l])
         FOR d1_i = 0l,d1_i_max-d1_i_min,1l DO BEGIN
           FOR d2_i = 0l,d2_i_max-d2_i_min,1l DO BEGIN
             IF (d3_sum EQ 1) THEN BEGIN
               sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d3_i_min:d3_i_max,d1_i_min+d1_i,d2_i_min+d2_i])
             ENDIF ELSE BEGIN
               sca_2dp[d1_i,d2_i] = sca_3dv[d3_i,d1_i_min+d1_i,d2_i_min+d2_i]
             ENDELSE
             
             d1_c = ds_ext.y[d1_i_min+d1_i] - ds_ext.y[d1_i_org]
             d2_c = ds_ext.z[d2_i_min+d2_i] - ds_ext.z[d2_i_org]      
             theta = (d1_c EQ 0 AND d2_c EQ 0) ? 0.0d0 : ATAN2(d2_c,d1_c)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360.0d0 - ABS(theta)
             ENDIF      
             theta_2dp[d1_i,d2_i] = theta
             rad_2dp[d1_i,d2_i] = SQRT(d1_c^2 + d2_c^2)    
           ENDFOR
         ENDFOR
       END
    1: BEGIN         
         d1_min = ds_ext.x[d1_i_min]
         d1_max = ds_ext.x[d1_i_max]
         d2_min = ds_ext.z[d2_i_min]
         d2_max = ds_ext.z[d2_i_max]
         dd1 = ABS(ds_ext.x[d1_i_min]-ds_ext.x[d1_i_min+1l])
         dd2 = ABS(ds_ext.z[d2_i_min]-ds_ext.z[d2_i_min+1l])
         FOR d1_i = 0l,d1_i_max-d1_i_min,1l DO BEGIN
           FOR d2_i = 0l,d2_i_max-d2_i_min,1l DO BEGIN
             IF (d3_sum EQ 1) THEN BEGIN
               sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d1_i_min+d1_i,d3_i_min:d3_i_max,d2_i_min+d2_i])
             ENDIF ELSE BEGIN
               sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i,d3_i,d2_i_min+d2_i]
             ENDELSE
             
             d1_c = ds_ext.x[d1_i_min+d1_i] - ds_ext.x[d1_i_org]
             d2_c = ds_ext.z[d2_i_min+d2_i] - ds_ext.z[d2_i_org]      
             theta = (d1_c EQ 0 AND d2_c EQ 0) ? 0.0d0 : ATAN2(d2_c,d1_c)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360.0d0 - ABS(theta)
             ENDIF      
             theta_2dp[d1_i,d2_i] = theta
             rad_2dp[d1_i,d2_i] = SQRT(d1_c^2 + d2_c^2)
           ENDFOR
         ENDFOR         
       END
    ELSE: BEGIN
            d1_min = ds_ext.x[d1_i_min]
            d1_max = ds_ext.x[d1_i_max]
            d2_min = ds_ext.y[d2_i_min]
            d2_max = ds_ext.y[d2_i_max]
                        
            dd1 = ABS(ds_ext.x[d1_i_min]-ds_ext.x[d1_i_min+1l])
            dd2 = ABS(ds_ext.y[d2_i_min]-ds_ext.y[d2_i_min+1l])
            FOR d1_i = 0l,d1_i_max-d1_i_min,1l DO BEGIN
              FOR d2_i = 0l,d2_i_max-d2_i_min,1l DO BEGIN
                IF (d3_sum EQ 1) THEN BEGIN
                  sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d1_i_min+d1_i,d2_i_min+d2_i,d3_i_min:d3_i_max])
                ENDIF ELSE BEGIN
                  sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i,d2_i_min+d2_i,d3_i]
                ENDELSE
                
                d1_c = ds_ext.x[d1_i_min+d1_i] - ds_ext.x[d1_i_org]
                d2_c = ds_ext.y[d2_i_min+d2_i] - ds_ext.y[d2_i_org]      
                theta = (d1_c EQ 0.0d0 AND d2_c EQ 0.0d0) ? 0.0d0 : ATAN2(d2_c,d1_c)
                theta = theta*180.0d0/!dpi
                IF (theta LT 0.0d0) THEN BEGIN
                  theta = 360.0d0 - ABS(theta)
                ENDIF      
                theta_2dp[d1_i,d2_i] = theta
                rad_2dp[d1_i,d2_i] = SQRT(d1_c^2 + d2_c^2)
              ENDFOR
            ENDFOR 
          END
  ENDCASE
  
  
  sel = WHERE(ABS(theta_2dp-theta_sel) LE theta_tol, sel_cnt)
  
  
  IF (sel_cnt GT 0l) THEN BEGIN
    IF (sca_min EQ sca_max) THEN BEGIN
      sca_min = MIN(sca_2dp[sel])
      sca_max = MAX(sca_2dp[sel])
    ENDIF
    
    IF STRLEN(output_filename) GT 0 THEN BEGIN
      PS_Start, output_filename, /ENCAPSULATED
    ENDIF 
    
    ;PLOT, rad_2dp[sel], sca_2dp[sel], psym=-2, xtitle=rad_title, ytitle=sca_title, xr=[0.0d0,rad_max], yr=[sca_min,sca_max]
    ;OPLOT, rad_2dp[sel], sca_2dp[sel], psym=-2
        
    IF STRLEN(output_filename) GT 0 THEN BEGIN  
      PS_End
    ENDIF
  ENDIF ELSE BEGIN
    PRINT, 'Error: no data.'
  ENDELSE
  
  IF STRLEN(data_filename) GT 0 THEN BEGIN    
    data_lun = 11
    OPENW, data_lun, data_filename
    FOR i = 0l,sel_cnt-1l,1l DO BEGIN
      PRINTF, data_lun, rad_2dp[sel[i]], ' ', sca_2dp[sel[i]]
    ENDFOR
    CLOSE, data_lun
  ENDIF
  
  
END
