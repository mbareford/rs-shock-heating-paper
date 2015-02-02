FUNCTION calc_vmag, crkdat, i
  RETURN, SQRT(crkdat.vx[i]^2+crkdat.vy[i]^2+crkdat.vz[i]^2)
END

FUNCTION calc_bmag, crkdat, i
  RETURN, SQRT(crkdat.bx[i]^2+crkdat.by[i]^2+crkdat.bz[i]^2)
END

FUNCTION calc_jmag, crkdat, i
  RETURN, SQRT(crkdat.jx[i]^2+crkdat.jy[i]^2+crkdat.jz[i]^2)
END

FUNCTION calc_fmag, crkdat, i
  fx = (crkdat.jy[i]*crkdat.bz[i] - crkdat.jz[i]*crkdat.by[i])
  fy = (crkdat.jz[i]*crkdat.bx[i] - crkdat.jx[i]*crkdat.bz[i])
  fz = (crkdat.jx[i]*crkdat.by[i] - crkdat.jy[i]*crkdat.bx[i])
          
  RETURN, SQRT(fx^2+fy^2+fz^2)
END

FUNCTION calc_emag, crkdat, i
  ex = crkdat.eta[i]*crkdat.jx[i] - (crkdat.vy[i]*crkdat.bz[i] - crkdat.vz[i]*crkdat.by[i])
  ey = crkdat.eta[i]*crkdat.jy[i] - (crkdat.vz[i]*crkdat.bx[i] - crkdat.vx[i]*crkdat.bz[i])
  ez = crkdat.eta[i]*crkdat.jz[i] - (crkdat.vx[i]*crkdat.by[i] - crkdat.vy[i]*crkdat.bx[i])
                
  RETURN, SQRT(ex^2+ey^2+ez^2)
END

FUNCTION calc_vxbmag, crkdat, i
  vxbx = (crkdat.vy[i]*crkdat.bz[i] - crkdat.vz[i]*crkdat.by[i])
  vxby = (crkdat.vz[i]*crkdat.bx[i] - crkdat.vx[i]*crkdat.bz[i])
  vxbz = (crkdat.vx[i]*crkdat.by[i] - crkdat.vy[i]*crkdat.bx[i])
  
  RETURN, SQRT(vxbx^2+vxby^2+vxbz^2)
END
          

PRO get_cork_plot_data, input_path, nx, ny, nz, xorg, yorg, zorg

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
    
  output_path = input_path +'/'

  lun0 = 10
  lun1 = 11
  lun2 = 12
  
  ; set up the cell centre arrays
  xc = DBLARR(nx+1)
  yc = DBLARR(ny+1)
  zc = DBLARR(nz+1)  
  dx = 2.0d0/DOUBLE(nx)
  dy = 2.0d0/DOUBLE(ny)
  dz = 20.0d0/DOUBLE(nz) 
  xc[0] = -2.0d0 - dx/2.0
  FOR i=1,nx,1 DO BEGIN
    xc[i] = -2.0d0 + dx*DOUBLE(i-1) + dx/2.0
  ENDFOR
  yc[0] = -2.0d0 - dy/2.0
  FOR i=1,ny,1 DO BEGIN
    yc[i] = -2.0d0 + dy*DOUBLE(i-1) + dy/2.0
  ENDFOR
  zc[0] = -20.0d0 - dz/2.0  
  FOR i=1,nz,1 DO BEGIN
    zc[i] = -20.0d0 + dz*DOUBLE(i-1) + dz/2.0
  ENDFOR
     
  crkdat = getcorkdata(wkdir=input_path)
  

  ; threshold set at background temperature
  ;visc_thres = 0.0d0
  visc_thres = 6.506d-5
      
    
  OPENW, lun0, output_path+'cork_ids.txt', WIDTH=1024
  OPENW, lun2, output_path+'cork_spds.txt', WIDTH=1024
  crkid = crkdat.id[UNIQ(crkdat.id, SORT(crkdat.id))]
  imax = SIZE(crkid,/N_ELEMENTS) - 1
  FOR i=LONG(0),imax,1 DO BEGIN
    
    PRINT, "Processing cork ", i+1, "..."
          
    xi = LONG(crkid[i]) MOD nx
    yi = (LONG(crkid[i])/nx) MOD ny
    zi = (LONG(crkid[i])/nx)/ny
    
    dx = ABS(xorg-xc[xi])
    dx = xorg GT xc[xi] ? -dx : dx
    dy = ABS(yorg-yc[yi])
    dy = yorg GT yc[yi] ? -dy : dy
    dz = ABS(zorg-zc[zi])
    dz = zorg GT zc[zi] ? -dz : dz
    
    rad = SQRT(dx^2 + dy^2 + dz^2)
    theta = rad GT 0.0d0 ? ACOS(dz/rad)*(180.0d0/!dpi) : 0.0d0
    phi = dx NE 0.0d0 ? ATAN(ABS(dy)/ABS(dx))*(180.0d0/!dpi) : 0.0d0
    IF (dy GT 0.0d0 && dx LT 0.0d0) THEN BEGIN
      phi = 90.0d0 + phi
    ENDIF ELSE BEGIN
      IF (dy LT 0.0d0 && dx LT 0.0d0) THEN BEGIN
        phi = 180.0d0 + phi
      ENDIF ELSE BEGIN
        IF (dy LT 0.0d0 && dx GT 0.0d0) THEN BEGIN
          phi = 270.0d0 + phi
        ENDIF
      ENDELSE
    ENDELSE
  
    ; write the starting positions of all the corks in the fleet
    PRINTF, lun0, LONG(crkid[i]), ' ', xc[xi], ' ', yc[yi], ' ', zc[zi], ' ', dx, ' ', dy, ' ', dz,  ' ', rad, ' ', theta, ' ', phi LT 0.0d0 ? phi + 360.0d0 : phi
   
    ; write data concerning a specific cork
    OPENW, lun1, output_path+'cork_'+STRTRIM(LONG(crkid[i]),2)+'.txt', WIDTH=1024
    
    ; find all the cork data entries for a specific cork
    cnt = 0l    
    sel = WHERE(crkid[i] EQ crkdat.id, cnt)
    
    IF (cnt GT 0l) THEN BEGIN
    
      imax2 = SIZE(sel,/N_ELEMENTS) - 1
      
      post_shock_tsc = -1l
      post_shock_t = -1.0d0
      shock_encountered = 0

      FOR i2=LONG(0),imax2,1 DO BEGIN
        ; iterate through the cork trail
        
        i3 = sel[i2]                                         
        i3m1 = i2 GT 0 ? sel[i2-1] : -1
        i3p1 = i2 LT imax2 ? sel[i2+1] : -1
          
        vmag = calc_vmag(crkdat, i3)
        bmag = calc_bmag(crkdat, i3)
        jmag = calc_jmag(crkdat, i3)
        fmag = calc_fmag(crkdat, i3)         
        emag = calc_emag(crkdat, i3)
        vxbmag = calc_vxbmag(crkdat, i3)
        
        ; determine the sound, alfven and magnetosonic speeds at a specific cork position        
        cs = SQRT((gamma-1.0d0)*crkdat.energy[i3]) 
        vA = bmag/SQRT(crkdat.rho[i3])
        cms2 = cs^2 + vA^2
        cms = SQRT(cms2)
        cms4 = cms2^2
        
        ; determine the slow and fast mode speeds at a specific cork position  
        vslow = 0.0d0
        vfast = 0.0d0
        vdotB = crkdat.vx[i3]*crkdat.bx[i3] + crkdat.vy[i3]*crkdat.by[i3] + crkdat.vz[i3]*crkdat.bz[i3]
        IF (vmag GT 0.0d0 AND bmag GT 0.0d0) THEN BEGIN
          aux = 4.0d0*(cs^2)*(vA^2)*((vdotB/(vmag*bmag))^2)          
          IF (cms4 GE aux) THEN BEGIN
            vfast = SQRT(0.5d0*(cms2 + SQRT(cms4 - aux)))            
            IF (cms2 GE SQRT(cms4 - aux)) THEN BEGIN
              vslow = SQRT(0.5d0*(cms2 - SQRT(cms4 - aux)))
            ENDIF
          ENDIF            
        ENDIF
                                                                                                
        ; calculate electric field              
        ex = crkdat.eta[i3]*crkdat.jx[i3] - (crkdat.vy[i3]*crkdat.bz[i3] - crkdat.vz[i3]*crkdat.by[i3])
        ey = crkdat.eta[i3]*crkdat.jy[i3] - (crkdat.vz[i3]*crkdat.bx[i3] - crkdat.vx[i3]*crkdat.bz[i3])
        ez = crkdat.eta[i3]*crkdat.jz[i3] - (crkdat.vx[i3]*crkdat.by[i3] - crkdat.vy[i3]*crkdat.bx[i3])
                                                                                                    
        PRINTF, lun1, LONG(crkdat.id[i3]), ' ', LONG(crkdat.new[i3]), ' ', crkdat.x[i3], ' ', crkdat.y[i3], ' ', crkdat.z[i3], ' ', $
                           crkdat.dt[i3], ' ', crkdat.t[i3], ' ', crkdat.ds[i3], ' ', crkdat.s[i3], ' ', $
                           crkdat.rho[i3], ' ', crkdat.energy[i3], ' ', crkdat.visc[i3], ' ', crkdat.ohmic[i3], ' ', $ 
                           crkdat.vx[i3], ' ', crkdat.vy[i3], ' ', crkdat.vz[i3], ' ', $
                           crkdat.bx[i3], ' ', crkdat.by[i3], ' ', crkdat.bz[i3], ' ', $
                           crkdat.jx[i3], ' ', crkdat.jy[i3], ' ', crkdat.jz[i3], ' ', $
                           ex, ' ', ey, ' ', ez, ' ', $
                           cs, ' ', vA, ' ', cms, ' ', vslow, ' ', vfast
        
        IF (vmag LT vslow) THEN BEGIN
          IF (shock_encountered EQ 1) THEN BEGIN
            post_shock_tsc = post_shock_tsc + 1l
            post_shock_t = post_shock_t + crkdat.dt[i3]
          ENDIF
        ENDIF ELSE BEGIN
          post_shock_tsc = 0l
          post_shock_t = 0.0d0
          shock_encountered = 1
        ENDELSE

        PRINTF, lun2, LONG(crkdat.id[i3]), ' ', crkdat.t[i3], ' ', crkdat.dt[i3], ' ', $
                      crkdat.rho[i3], ' ', crkdat.energy[i3], ' ', crkdat.visc[i3], ' ', $
                      vmag, ' ', cs, ' ', vA, ' ', vslow, ' ', vfast, ' ', post_shock_tsc, ' ', post_shock_t  
      ENDFOR
      
    ENDIF
    CLOSE, lun1
      
  ENDFOR
  CLOSE, lun0
  CLOSE, lun2
      
END
