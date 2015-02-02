FUNCTION get_data_3dv, t_i,input_path,cen, data_id
  data_3dv = 0 
  
  mu0 = 4.0d-7*!dpi
  mp = 1.6726d-27
  me = 9.1094d-31
  ep = 8.8542d-12
  ec = 1.6022d-19
  kB = 1.3807d-23
  gamma = 5.0d0/3.0d0  
  mf = 1.0d0
  m0 = mf*mp

  B0 = 5.0d-3  ; T
  L0 = 1.0d6   ; m
  n0 = 1.0d15  ; m^-3
  rho0 = m0*n0            ; kg m^-3
  vA = B0/SQRT(mu0*rho0)  ; m s^-1
  tA = L0/vA              ; s

  P0 = B0^2/mu0     ; N m^-2
  En0 = P0/rho0     ; J 
  Tp0 = En0*(m0/kB) ; K  
  j0 = B0/(mu0*L0)  ; K
    
  
  CASE data_id OF    
    1: BEGIN
         ;;; density ;;;
         ds = getdata(t_i,wkdir=input_path,/rho)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho)           
         ENDELSE
         data_3dv = ds_ext.rho
       END 
    2: BEGIN
         ;;; density difference ;;;
         ds0 = getdata(0,wkdir=input_path,/rho)
         ds = getdata(t_i,wkdir=input_path,/rho)
         IF (cen EQ 0) THEN BEGIN
           ds0_ext = get_ext_ccdata(ds0,/rho)
           ds_ext = get_ext_ccdata(ds,/rho)
         ENDIF ELSE BEGIN
           ds0_ext = get_ext_vcdata(ds0,/rho)
           ds_ext = get_ext_vcdata(ds,/rho)           
         ENDELSE  
         data_3dv = ds_ext.rho-ds0_ext.rho                  
       END       
    3: BEGIN
         ;;; temperature ;;;
         ds = getdata(t_i,wkdir=input_path,/energy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/temperature)           
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/temperature)           
         ENDELSE  
         data_3dv = ds_ext.temperature
         data_3dv = ds_ext.temperature*(Tp0/1.0d6)
       END             
    4: BEGIN
         ;;; temperature difference ;;;
         ds0 = getdata(0,wkdir=input_path,/temperature,/energy)
         ds = getdata(t_i,wkdir=input_path,/temperature,/energy)
         IF (cen EQ 0) THEN BEGIN
           ds0_ext = get_ext_ccdata(ds0,/temperature)
           ds_ext = get_ext_ccdata(ds,/temperature)
         ENDIF ELSE BEGIN
           ds0_ext = get_ext_vcdata(ds0,/temperature)
           ds_ext = get_ext_vcdata(ds,/temperature)           
         ENDELSE  
         data_3dv = (ds_ext.temperature-ds0_ext.temperature)*(Tp0/1.0d6)
       END
    5: BEGIN
         ;;; eta ;;;
         ds = getdata(t_i,wkdir=input_path,/eta)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/resistivity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/resistivity)           
         ENDELSE  
         data_3dv = ds_ext.eta
       END
    6: BEGIN
         ;;; energy ;;;
         ds = getdata(t_i,wkdir=input_path,/energy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/energy)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/energy)           
         ENDELSE  
         data_3dv = (ds_ext.energy/3.0d0)*(Tp0/1.0d6)
         sel_cnt = 0l
         sel_i = WHERE(data_3dv LE 0.02d0, sel_cnt)
         IF (sel_cnt GT 0L) THEN BEGIN
           data_3dv[sel_i] = 0.0d0
         ENDIF
       END             
    7: BEGIN
         ;;; energy difference ;;;
         ds0 = getdata(0,wkdir=input_path,/energy)
         ds = getdata(t_i,wkdir=input_path,/energy)
         IF (cen EQ 0) THEN BEGIN
           ds0_ext = get_ext_ccdata(ds0,/energy)
           ds_ext = get_ext_ccdata(ds,/energy)
         ENDIF ELSE BEGIN
           ds0_ext = get_ext_vcdata(ds0,/energy)
           ds_ext = get_ext_vcdata(ds,/energy)           
         ENDELSE  
         data_3dv = ((ds_ext.energy-ds0_ext.energy)/3.0d0)*(Tp0/1.0d6)
       END 
       
    8: BEGIN
         ;;; jcrit ;;;
         jcrit_const = 75840.0d0
         ds = getdata(t_i,wkdir=input_path,/rho,/energy,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/temperature,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/temperature,/bfield)           
         ENDELSE  
         data_3dv = jcrit_const*ds_ext.rho*ds_ext.temperature/SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
       END 

    9: BEGIN
         ;;; jmag ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/current)           
         ENDELSE  
         data_3dv = SQRT(ds_ext.jx^2 + ds_ext.jy^2 + ds_ext.jz^2)         
       END
   
   10: BEGIN
         ;;; jsuper ;;;
         jcrit_const = 75840.0d0
         ds = getdata(t_i,wkdir=input_path,/rho,/energy,/bx,/by,/bz)
         nx = ds.grid.npts[0]
         ny = ds.grid.npts[1]
         nz = ds.grid.npts[2]
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/temperature,/bfield,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/temperature,/bfield,/current)           
         ENDELSE
         data_3dv = SQRT(ds_ext.jx^2 + ds_ext.jy^2 + ds_ext.jz^2)
         jcrit = jcrit_const*ds_ext.rho*ds_ext.temperature/SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
         
         cnt = 0l
         sel = WHERE(data_3dv GE jcrit, cnt)
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = 1.0d0
           
           x_i = sel mod nx
           y_i = (sel/nx) mod ny
           z_i = sel/(nx*ny)
           
           FOR i = 0,cnt-1,1 DO BEGIN
             PRINT, "x=", ds.grid.x[x_i[i]], ", y=", ds.grid.y[y_i[i]], ", z=", ds.grid.z[z_i[i]], "."
           ENDFOR
         ENDIF
         PRINT, "There are ", cnt, " cells containing supercritical current."
         sel = WHERE(data_3dv LT jcrit, cnt)
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = 0.0d0
         ENDIF
       END
        
   11: BEGIN
         ;;; ohmic_heat ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/ohmic_heat)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/ohmic_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/ohmic_heat)
         ENDELSE
         data_3dv = ds_ext.ohmic_heat/ds_ext.rho
       END              
   12: BEGIN
         ;;; visc_heat ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat)        
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/visc_heat)           
         ENDELSE
         data_3dv = ds_ext.visc_heat/ds_ext.rho
       END
   13: BEGIN
         ;;; visc_heat_xy ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xy_visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xy_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.xy_visc_heat/ds_ext.rho
       END
   14: BEGIN
         ;;; visc_heat_xz ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xz_visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xz_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.xz_visc_heat/ds_ext.rho
       END
   15: BEGIN
         ;;; visc_heat_yz ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_yz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/yz_visc_heat)           
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/yz_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.yz_visc_heat/ds_ext.rho
       END
   16: BEGIN
         ;;; visc_heat_xx ;;; 
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_xx)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/xx_visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/xx_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.xx_visc_heat/ds_ext.rho
       END
   17: BEGIN
         ;;; visc_heat_yy ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_yy)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/yy_visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/yy_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.yy_visc_heat/ds_ext.rho
       END
   18: BEGIN
         ;;; visc_heat_zz ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/visc_heat_zz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/zz_visc_heat)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/zz_visc_heat)           
         ENDELSE
         data_3dv = ds_ext.zz_visc_heat/ds_ext.rho
       END
       
     
   21: BEGIN
         ;;; |B| ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)
         ENDELSE
         data_3dv = SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
       END    
   22: BEGIN
         ;;; |v| ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity)
         ENDELSE
         data_3dv = SQRT(ds_ext.vx^2 + ds_ext.vy^2 + ds_ext.vz^2)
       END
   23: BEGIN
         ;;; |j| ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current)
         ENDELSE
         data_3dv = SQRT(ds_ext.jx^2 + ds_ext.jy^2 + ds_ext.jz^2)                  
       END    
   24: BEGIN
         ;;; |f| ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity,/bfield,/current,/lorentz)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity,/bfield,/current,/lorentz)
         ENDELSE
         data_3dv = SQRT(ds_ext.fx^2 + ds_ext.fy^2 + ds_ext.fz^2)
       END
   25: BEGIN
         ;;; |lagrangian f| ;;;
         ds = getdata(t_i,wkdir=input_path,/fx,/fy,/fz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/ffield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/ffield)
         ENDELSE
         data_3dv = SQRT(ds_ext.fx^2 + ds_ext.fy^2 + ds_ext.fz^2)
         
         cnt = 0l
         sel = WHERE(data_3dv LT 0.01d0, cnt)
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = 0.0d0
         ENDIF
       END
   26: BEGIN
         ;;; |vorticity| ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/vorticity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/vorticity)
         ENDELSE
         data_3dv = SQRT(ds_ext.wx^2 + ds_ext.wy^2 + ds_ext.wz^2)
       END
   27: BEGIN
         ;;; x component of vorticity ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/vorticity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/vorticity)
         ENDELSE
         data_3dv = ds_ext.wx
       END
   28: BEGIN
         ;;; y component of vorticity ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/vorticity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/vorticity)
         ENDELSE
         data_3dv = ds_ext.wy
       END
   29: BEGIN
         ;;; z component of vorticity ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/vorticity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/vorticity)
         ENDELSE
         data_3dv = ds_ext.wz
       END
          
                 
   31: BEGIN
         ;;; bz ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)           
         ENDELSE
         data_3dv = ds_ext.bz
       END
   32: BEGIN
         ;;; vz ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity)           
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity)
         ENDELSE
         data_3dv = ds_ext.vz
       END
   33: BEGIN
         ;;; jz ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current)           
         ENDELSE
         data_3dv = ds_ext.jz
       END
   34: BEGIN
         ;;; fz ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current,/lorentz)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current,/lorentz)           
         ENDELSE
         data_3dv = ds_ext.fz
       END   
   35: BEGIN
         ;;; 'lagrangian fz' ;;;
         ds = getdata(t_i,wkdir=input_path,/fx,/fy,/fz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/ffield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/ffield)           
         ENDELSE
         data_3dv = ds_ext.fz
       END                               
   36: BEGIN
         ;;; btheta ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)           
         ENDELSE
         
         nx = SIZE(ds.grid.x,/N_ELEMENTS) 
         ny = SIZE(ds.grid.y,/N_ELEMENTS)
         data_3dv = ds_ext.bz
  
         FOR x_i = 0,nx-1l,1l DO BEGIN
           FOR y_i = 0,ny-1l,1l DO BEGIN
             x = ds_ext.x[x_i]
             y = ds_ext.y[y_i]      
         
             theta = (x EQ 0 AND y EQ 0) ? 0.0d0 : ATAN2(y,x)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360.0d0 - ABS(theta)
             ENDIF
             theta = theta*!dpi/180.0d0
         
             bx = ds_ext.bx[x_i,y_i,*]
             by = ds_ext.by[x_i,y_i,*]
             data_3dv[x_i,y_i,*] = ds_ext.by[x_i,y_i,*]*COS(theta) - ds_ext.bx[x_i,y_i,*]*SIN(theta)
           ENDFOR        
         ENDFOR         
       END 
   37: BEGIN
         ;;; br ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)           
         ENDELSE
         
         nx = SIZE(ds.grid.x,/N_ELEMENTS) 
         ny = SIZE(ds.grid.y,/N_ELEMENTS)
         data_3dv = ds_ext.bz
  
         FOR x_i = 0,nx-1l,1l DO BEGIN
           FOR y_i = 0,ny-1l,1l DO BEGIN
             x = ds_ext.x[x_i]
             y = ds_ext.y[y_i]      
         
             theta = (x EQ 0 AND y EQ 0) ? 0.0d0 : ATAN2(y,x)
             theta = theta*180.0d0/!dpi
             IF (theta LT 0.0d0) THEN BEGIN
               theta = 360.0d0 - ABS(theta)
             ENDIF
             theta = theta*!dpi/180.0d0
         
             bx = ds_ext.bx[x_i,y_i,*]
             by = ds_ext.by[x_i,y_i,*]
             data_3dv[x_i,y_i,*] = ds_ext.by[x_i,y_i,*]*COS(theta) + ds_ext.bx[x_i,y_i,*]*SIN(theta)
           ENDFOR        
         ENDFOR         
       END
   38: BEGIN
         ;;; vx ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity)           
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity)
         ENDELSE
         data_3dv = ds_ext.vx
       END             
                                          
   41: BEGIN
         ;;; j.B ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity,/bfield,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity,/bfield,/current)
         ENDELSE
         bmag = SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
         data_3dv = ABS(ds_ext.jx*ds_ext.bx + ds_ext.jy*ds_ext.by + ds_ext.jz*ds_ext.bz)/bmag
       END      
       
         
   42: BEGIN
         ;;; (j.B)/(|j||B|) ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity,/bfield)
         ENDELSE
         bmag = SQRT(ds_ext.bx^2+ds_ext.by^2+ds_ext.bz^2)
         vmag = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)
         data_3dv = ABS(ds_ext.vx*ds_ext.bx + ds_ext.vy*ds_ext.by + ds_ext.vz*ds_ext.bz)/(vmag*bmag)
       END              
   43: BEGIN
         ;;; e.B ;;;
         ds = getdata(t_i,wkdir=input_path,/eta,/bx,/by,/bz,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/efield,/bfield,/current,/resistivity,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/efield,/bfield,/current,/resistivity,/velocity)
         ENDELSE
         data_3dv = ds_ext.ex*ds_ext.bx + ds_ext.ey*ds_ext.by + ds_ext.ez*ds_ext.bz         
       END
   
   
   51: BEGIN
         ;;; sound speed mach number ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/pressure,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/pressure,/velocity)
         ENDELSE
         data_3dv = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)/SQRT(5.0d0/3.0d0*ds_ext.p/ds_ext.rho)
         nonshock_cnt = 0l
         nonshock_i = WHERE(data_3dv LT 1.0d0, nonshock_cnt)
         IF (nonshock_cnt GT 0l) THEN BEGIN
           data_3dv[nonshock_i] = 0.0d0
         ENDIF
       END           
   52: BEGIN
         ;;; vA mach number ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/bfield)
         ENDELSE
         data_3dv = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)/(SQRT(ds_ext.bx^2+ds_ext.by^2+ds_ext.bz^2)/SQRT(ds_ext.rho))         
         nonshock_cnt = 0l
         nonshock_i = WHERE(data_3dv LT 1.0d0, nonshock_cnt)
         IF (nonshock_cnt GT 0l) THEN BEGIN
           data_3dv[nonshock_i] = 0.0d0
         ENDIF
       END
   53: BEGIN
         ;;; cms ;;; 
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/pressure,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/pressure,/bfield)
         ENDELSE
         data_3dv = SQRT(((5.0d0/3.0d0)*(ds_ext.p/ds_ext.rho)) + ((ds_ext.bx^2+ds_ext.by^2+ds_ext.bz^2)/ds_ext.rho))
       END
          
   54: BEGIN
         ;;; slow speed shock ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature,/visc_heat,/bx,/by,/bz,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/pressure,/visc_heat,/bfield,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/pressure,/visc_heat,/bfield,/velocity)
         ENDELSE
         
         bmag = SQRT(ds_ext.bx^2+ds_ext.by^2+ds_ext.bz^2)
         vmag = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)                  
         
         cs = SQRT(((5.0d0/3.0d0)*ds_ext.p)/ds_ext.rho)
         va = bmag/SQRT(ds_ext.rho)
         cms2 = cs^2 + va^2
         cms = SQRT(cms2)
         cms4 = cms2*cms2
         
         vdotB = ds_ext.vx*ds_ext.bx+ds_ext.vy*ds_ext.by+ds_ext.vz*ds_ext.bz
         aux = 4.0d0*(cs^2)*(va^2)*(vdotB^2)
         
         vslow = cs
         
         cnt = 0l
         sel = WHERE(vmag GT 0.0d0, cnt)
         IF (cnt GT 0l) THEN BEGIN        

           aux[sel] = aux[sel]/((vmag[sel]^2)*(bmag[sel]^2))
           
           cnt2 = 0l
           sel2 = WHERE(cms4[sel] GT aux[sel], cnt2)
           IF (cnt2 GT 0l) THEN BEGIN
             sel2i = sel[sel2]             
             cnt3 = 0l
             sel3 = WHERE(cms2[sel2i] GT SQRT(cms4[sel2i] - aux[sel2i]), cnt3)
             IF (cnt3 GT 0l) THEN BEGIN
               sel3i = sel2i[sel3]
               vslow[sel3i] = SQRT(0.5d0*(cms2[sel3i] - SQRT(cms4[sel3i] - aux[sel3i])))               
             ENDIF
           ENDIF
           
           cnt2 = 0l
           sel2 = WHERE(cms4[sel] LE aux[sel], cnt2)
           IF (cnt2 GT 0l) THEN BEGIN
             vmag[sel[sel2]] = 0.0d0
           ENDIF
           
           cnt2 = 0l
           sel2 = WHERE(cms2[sel] LE SQRT(cms4[sel] - aux[sel]), cnt2)
           IF (cnt2 GT 0l) THEN BEGIN
             vmag[sel[sel2]] = 0.0d0
           ENDIF
         ENDIF      
             
         data_3dv = vmag/vslow
                   
         cnt = 0l
         sel = WHERE(data_3dv LT 1.0d0, cnt)
         PRINT,'sub slow cnt=',cnt
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = 0.0d0
         ENDIF                          
         
         cnt = 0l
         sel = WHERE(data_3dv GT 0.0d0, cnt)
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = ALOG10(data_3dv[sel])
         ENDIF 
                                 
       END        
   55: BEGIN
         ;;; fast speed shock ;;;
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature,/visc_heat,/bx,/by,/bz,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/rho,/pressure,/visc_heat,/bfield,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/rho,/pressure,/visc_heat,/bfield,/velocity)
         ENDELSE
         
         bmag = SQRT(ds_ext.bx^2+ds_ext.by^2+ds_ext.bz^2)
         vmag = SQRT(ds_ext.vx^2+ds_ext.vy^2+ds_ext.vz^2)
         
         cs = SQRT(((5.0d0/3.0d0)*ds_ext.p)/ds_ext.rho)
         va = bmag/SQRT(ds_ext.rho)
         cms2 = cs^2 + va^2
         cms = SQRT(cms2)
         cms4 = cms2*cms2
         
         vdotB = ds_ext.vx*ds_ext.bx+ds_ext.vy*ds_ext.by+ds_ext.vz*ds_ext.bz
         aux = 4.0d0*(cs^2)*(va^2)*(vdotB^2)
         
         vfast = va
         
         cnt = 0l
         sel = WHERE(vmag GT 0.0d0, cnt)
         IF (cnt GT 0l) THEN BEGIN        

           aux[sel] = aux[sel]/((vmag[sel]^2)*(bmag[sel]^2))
           
           sel2 = WHERE(cms4[sel] GT aux[sel], cnt2)
           IF (cnt2 GT 0l) THEN BEGIN
             sel2i = sel[sel2]             
             vfast[sel2i] = SQRT(0.5d0*(cms2[sel2i] + SQRT(cms4[sel2i] - aux[sel2i])))                          
           ENDIF
           
           sel2 = WHERE(cms4[sel] LE aux[sel], cnt2)
           IF (cnt2 GT 0l) THEN BEGIN
             vmag[sel[sel2]] = 0.0d0
           ENDIF                      
         ENDIF                        
         
         data_3dv = vmag/vfast
                  
         cnt = 0l
         sel = WHERE(data_3dv LT 1.0d0, cnt)
         PRINT,'sub fast cnt=',cnt
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = 0.0d0
         ENDIF                          
         
         cnt = 0l
         sel = WHERE(data_3dv GT 0.0d0, cnt)
         IF (cnt GT 0l) THEN BEGIN
           data_3dv[sel] = ALOG10(data_3dv[sel])
         ENDIF                  
       END 
        
    61: BEGIN
         ;;; pressure force ;;;
         ds = getdata(t_i,wkdir=input_path,/temperature,/rho)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/pressure)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/pressure)           
         ENDELSE  
         data_3dv = -SQRT(ds_ext.dxp^2 + ds_ext.dyp^2 + ds_ext.dzp^2)
       END                  
    
    62: BEGIN
         ;;; |f| ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity,/bfield,/current,/lorentz)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity,/bfield,/current,/lorentz)
         ENDELSE
         data_3dv = SQRT(ds_ext.fx^2 + ds_ext.fy^2 + ds_ext.fz^2)/SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
       END
    
    63: BEGIN
         ;;; pressure ;;;
         ds = getdata(t_i,wkdir=input_path,/temperature,/rho)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/pressure)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/pressure)           
         ENDELSE  
         data_3dv = ds_ext.p
       END
    
    64: BEGIN
         ;;; pressure difference ;;;
         ds0 = getdata(0,wkdir=input_path,/rho,/temperature)
         ds = getdata(t_i,wkdir=input_path,/rho,/temperature)
         IF (cen EQ 0) THEN BEGIN
           ds0_ext = get_ext_ccdata(ds0,/pressure)
           ds_ext = get_ext_ccdata(ds,/pressure)
         ENDIF ELSE BEGIN
           ds0_ext = get_ext_vcdata(ds0,/pressure)
           ds_ext = get_ext_vcdata(ds,/pressure)           
         ENDELSE  
         data_3dv = ds_ext.p-ds0_ext.p
       END  
       
    65: BEGIN
         ;;; jmag difference ;;;
         ds0 = getdata(0,wkdir=input_path,/bx,/by,/bz)
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds0_ext = get_ext_ccdata(ds0,/current)
           ds_ext = get_ext_ccdata(ds,/current)
         ENDIF ELSE BEGIN
           ds0_ext = get_ext_vcdata(ds0,/current)
           ds_ext = get_ext_vcdata(ds,/current)
         ENDELSE
         jmag0 = SQRT(ds0_ext.jx^2 + ds0_ext.jy^2 + ds0_ext.jz^2)
         jmag = SQRT(ds_ext.jx^2 + ds_ext.jy^2 + ds_ext.jz^2)
         data_3dv = jmag - jmag0
       END

    66: BEGIN
         ;;; plasma beta ;;;
         ds = getdata(t_i,wkdir=input_path,/temperature,/rho,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/pressure,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/pressure,/bfield)           
         ENDELSE  
         data_3dv = 2.0d0*ds_ext.p / SQRT(ds_ext.bx^2 + ds_ext.by^2 + ds_ext.bz^2)
       END
                                                         
    ELSE: BEGIN
           PRINT, 'Error: unrecognised data id.'
           RETURN, data_3dv           
         END
  ENDCASE
  
  RETURN, data_3dv
END

