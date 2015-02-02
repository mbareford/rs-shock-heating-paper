; Return an extended data set for a particular snapshot.
;
; The data set comprises the following.
;
; Cell boundary coordinates.
;
; Cell-centred pressures and face-centred pressure gradients.
;
; Face-centred magnetic field components.
; Vertex-centred velocity components.
; Edge-centred currents, resistivities and vorticities (calculated from face-centred velocities).
;
; Electric field components calculated from edge-centred currents and resistivities,
; face-centred magnetic field components, and vertex-centred velocity components.
;
; Alpha values calculated from edge-centred currents and face-centred magnetic field components.
;
; ds: the basic data set for some snapshot
; jcrit: the current threshold above which anomalous resistivity is applied
; eta0: the background resistivity
; eta1: the anomalous resistivity
FUNCTION get_ext_data, ds, jcrit, eta0, eta1

  ON_ERROR, 2
  CLOSE, 1

  IF (N_ELEMENTS(ds) NE 0) THEN BEGIN
    
    prec = ds.grid.prec
    unit = (prec EQ 4) ? 1.0 : 1.0D
    zero = (prec EQ 4) ? 0.0 : 0.0D

    
    nx = ds.grid.npts[0]
    ny = ds.grid.npts[1]
    nz = ds.grid.npts[2]
   
    ; calculate cell dimensions
    dx = SHIFT(ds.grid.x,-1) - ds.grid.x
    dy = SHIFT(ds.grid.y,-1) - ds.grid.y
    dz = SHIFT(ds.grid.z,-1) - ds.grid.z
    
    ; correct the last element in di (i=x,y,z) arrays
    dx[nx-1] = ds.grid.x[nx-1] - ds.grid.x[nx-2]
    dy[ny-1] = ds.grid.y[ny-1] - ds.grid.y[ny-2]
    dz[nz-1] = ds.grid.z[nz-1] - ds.grid.z[nz-2]
    
    ; convert di to 3D arrays                   
    adx = REFORM((dx # REPLICATE(unit,ny))[*] # REPLICATE(unit,nz), nx,ny,nz)
    ady = REFORM((REPLICATE(unit,nx) # dy)[*] # REPLICATE(unit,nz), nx,ny,nz)
    adz = REFORM((REPLICATE(unit,nx) # REPLICATE(unit,ny))[*] # dz, nx,ny,nz)


    ;print, 'using cell boundary coordinates...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   


    ;print, 'using vertex-centred velocities...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    ;print, 'calculating cell-centred pressure and face-centred pressure gradients...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; calculate the thermal pressure
    p = ds.rho*ds.temperature

    ; calculate the gradient of the thermal pressure
    dxp = (SHIFT(p,-1,0,0) - p) / adx
    dyp = (SHIFT(p,0,-1,0) - p) / ady
    dzp = (SHIFT(p,0,0,-1) - p) / adz

    ; enforce zero pressure gradient at maximum boundary    
    dxp[nx-1,*,*] = dxp[nx-2,*,*]
    dxp[*,ny-1,*] = dxp[*,ny-2,*]
    dxp[*,*,nz-1] = dxp[*,*,nz-2]
    dyp[nx-1,*,*] = dyp[nx-2,*,*]
    dyp[*,ny-1,*] = dyp[*,ny-2,*]
    dyp[*,*,nz-1] = dyp[*,*,nz-2]
    dzp[nx-1,*,*] = dzp[nx-2,*,*]
    dzp[*,ny-1,*] = dzp[*,ny-2,*]
    dzp[*,*,nz-1] = dzp[*,*,nz-2]     
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    ;print, 'calculating face-centred velocities...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    vxf = ds.vx
    vyf = ds.vy
    vzf = ds.vz

    vxf[1:nx-1,1:ny-1,1:nz-1] = (ds.vx[1:nx-1,1:ny-1,1:nz-1] + ds.vx[1:nx-1,0:ny-2,1:nz-1] + ds.vx[1:nx-1,1:ny-1,0:nz-2] + ds.vx[1:nx-1,0:ny-2,0:nz-2]) / 4.0
    vyf[1:nx-1,1:ny-1,1:nz-1] = (ds.vy[1:nx-1,1:ny-1,1:nz-1] + ds.vy[0:nx-2,1:ny-1,1:nz-1] + ds.vy[1:nx-1,1:ny-1,0:nz-2] + ds.vy[0:nx-2,1:ny-1,0:nz-2]) / 4.0
    vzf[1:nx-1,1:ny-1,1:nz-1] = (ds.vz[1:nx-1,1:ny-1,1:nz-1] + ds.vz[0:nx-2,1:ny-1,1:nz-1] + ds.vz[1:nx-1,0:ny-2,1:nz-1] + ds.vz[0:nx-2,0:ny-2,1:nz-1]) / 4.0

    vxf[1:nx-1,0,1:nz-1] = (ds.vx[1:nx-1,0,1:nz-1] + ds.vx[1:nx-1,0,0:nz-2]) / 2.0    
    vxf[1:nx-1,1:ny-1,0] = (ds.vx[1:nx-1,1:ny-1,0] + ds.vx[1:nx-1,0:ny-2,0]) / 2.0    
    vyf[0,1:ny-1,1:nz-1] = (ds.vy[0,1:ny-1,1:nz-1] + ds.vy[0,1:ny-1,0:nz-2]) / 2.0    
    vyf[1:nx-1,1:ny-1,0] = (ds.vy[1:nx-1,1:ny-1,0] + ds.vy[0:nx-2,1:ny-1,0]) / 2.0    
    vzf[0,1:ny-1,1:nz-1] = (ds.vz[0,1:ny-1,1:nz-1] + ds.vz[0,0:ny-2,1:nz-1]) / 2.0
    vzf[1:nx-1,0,1:nz-1] = (ds.vz[1:nx-1,0,1:nz-1] + ds.vz[0:nx-2,0,1:nz-1]) / 2.0
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
    ;print, 'calculating edge-centred vorticities (using face-centred velocity components)...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    wxe = ((SHIFT(vzf,0,-1,0) - vzf) / ady) - ((SHIFT(vyf,0,0,-1) - vyf) / adz)
    wye = ((SHIFT(vxf,0,0,-1) - vxf) / adz) - ((SHIFT(vzf,-1,0,0) - vzf) / adx)
    wze = ((SHIFT(vyf,-1,0,0) - vyf) / adx) - ((SHIFT(vxf,0,-1,0) - vxf) / ady)

    ; enforce zero vorticity since velocities are zero at boundary
    wxe[nx-1,*,*] = zero
    wye[nx-1,*,*] = zero
    wze[nx-1,*,*] = zero
    wxe[*,ny-1,*] = zero
    wye[*,ny-1,*] = zero
    wze[*,ny-1,*] = zero
    wxe[*,*,nz-1] = zero
    wye[*,*,nz-1] = zero
    wze[*,*,nz-1] = zero
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    


    ;print, 'calculating edge-centred currents (using face-centred magnetic field components)...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    jxe = ((SHIFT(ds.bz,0,-1,0) - ds.bz) / ady) - ((SHIFT(ds.by,0,0,-1) - ds.by) / adz)
    jye = ((SHIFT(ds.bx,0,0,-1) - ds.bx) / adz) - ((SHIFT(ds.bz,-1,0,0) - ds.bz) / adx)
    jze = ((SHIFT(ds.by,-1,0,0) - ds.by) / adx) - ((SHIFT(ds.bx,0,-1,0) - ds.bx) / ady)

    ; enforce zero current since magnetic field gradient is zero at boundaries
    jxe[nx-1,*,*] = zero
    jye[nx-1,*,*] = zero
    jze[nx-1,*,*] = zero
    jxe[*,ny-1,*] = zero
    jye[*,ny-1,*] = zero
    jze[*,ny-1,*] = zero
    jxe[*,*,nz-1] = zero
    jye[*,*,nz-1] = zero
    jze[*,*,nz-1] = zero
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
    
          

    ;print, 'calculating edge-centred resistivities...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; initialise the 3D resistivity array
    etae = jxe
    etae[*,*,*] = eta0
    etae_anom_i = WHERE(SQRT(jxe^2 + jye^2 + jze^2) GE ABS(jcrit), cnt)
    IF cnt GT 0 THEN BEGIN
      x_i = etae_anom_i mod nx
      y_i = (etae_anom_i/nx) mod ny
      z_i = etae_anom_i/(nx*ny)      
      etae[x_i,y_i,z_i] = eta1      
    ENDIF
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    


    ;print, 'calculating electric field components...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ex = etae*jxe - (ds.vy*ds.bz - ds.vz*ds.by)
    ey = etae*jye - (ds.vz*ds.bx - ds.vx*ds.bz)
    ez = etae*jze - (ds.vx*ds.by - ds.vy*ds.bx)
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   


   
    ;print, 'calculating alpha...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    alpha = jxe
    alpha[*,*,*] = 0.0
    bmag_nz_i = WHERE((ds.bx^2 + ds.by^2 + ds.bz^2) GT 0.0, cnt)
    IF cnt GT 0 THEN BEGIN
      x_i = bmag_nz_i mod nx
      y_i = (bmag_nz_i/nx) mod ny
      z_i = bmag_nz_i/(nx*ny)
      alpha[x_i,y_i,z_i] = (ds.bx[x_i,y_i,z_i]*jxe[x_i,y_i,z_i] + ds.by[x_i,y_i,z_i]*jye[x_i,y_i,z_i] + ds.bz[x_i,y_i,z_i]*jze[x_i,y_i,z_i])$
                           /(ds.bx[x_i,y_i,z_i]^2 + ds.by[x_i,y_i,z_i]^2 + ds.bz[x_i,y_i,z_i]^2)
    ENDIF
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    ;print, 'calculating j x b...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nff = SQRT((jye*ds.bz - jze*ds.by)^2 + (jze*ds.bx - jxe*ds.bz)^2 + (jxe*ds.by - jye*ds.bx)^2)
    ;print, 'finished.'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                         
    ext_data = CREATE_STRUCT('x',ds.grid.x, 'y',ds.grid.y, 'z',ds.grid.z,$
                             'p',p, 'dxp',dxp, 'dyp',dyp, 'dzp',dzp,$
                             'vx',ds.vx, 'vy',ds.vy, 'vz',ds.vz,$
                             'bx',ds.bx, 'by',ds.by, 'bz',ds.bz,$
                             'jx',jxe, 'jy',jye, 'jz',jze,$
                             'wx',wxe, 'wy',wye, 'wz',wze,$
                             'alpha',alpha, 'nff',nff, 'eta',etae,$
                             'ex',ex, 'ey',ey, 'ez',ez)        

  ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data structure>"
    ext_data = 0
  ENDELSE  

  RETURN, ext_data

END
