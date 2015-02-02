; Return an extended data set for a particular snapshot.
;
; The data set comprises the following.
;
; Cell centre coordinates.
; Cell-centred density, energy, temperature and pressure.
; Cell-centred Ohmic heating.
; Cell-centred viscous heating components.
; Cell-centred magnetic field components.
; Cell-centred velocity components.
; Cell-centred currents, resistivities and vorticities (calculated from face-centred velocities).
; Cell-centred Lorentz forces.
; Cell-centred raw Lorentz forces.
;
; Cell-centred electric field components calculated from cell-centred currents and resistivities,
; Cell-centred magnetic field components, and cell-centred velocity components.
;
; Alpha values calculated from cell-centred currents and cell-centred magnetic field components.
;
; ds: the basic data set for some snapshot
FUNCTION get_ext_ccdata, ds, rho=rho, energy=energy, temperature=temperature, pressure=pressure, ohmic_heat=ohmic_heat, $
                             visc_heat=visc_heat, xy_visc_heat=xy_visc_heat, xz_visc_heat=xz_visc_heat, yz_visc_heat=yz_visc_heat, $
                             xx_visc_heat=xx_visc_heat, yy_visc_heat=yy_visc_heat, zz_visc_heat=zz_visc_heat, $
                             velocity=velocity, vorticity=vorticity, $
                             bfield=bfield, efield=efield, resistivity=resistivity, current=current, $
                             lorentz=lorentz, ffield=ffield, alpha=alpha, dtW=dtW, all=all
  ON_ERROR, 2
  CLOSE, 1

  IF (N_ELEMENTS(ds) NE 0) THEN BEGIN
    
    gamma = 5.0d0/3.0d0
    
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

    ;print, 'calculating cell-centred cell coordinates...'
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    xc = ds.grid.x
    yc = ds.grid.y
    zc = ds.grid.z
    
    xc[1:nx-1] = (ds.grid.x[0:nx-2] + ds.grid.x[1:nx-1])/2.0d0
    yc[1:ny-1] = (ds.grid.y[0:ny-2] + ds.grid.y[1:ny-1])/2.0d0
    zc[1:nz-1] = (ds.grid.z[0:nz-2] + ds.grid.z[1:nz-1])/2.0d0
    
    ; correct the first cell centre coordinates
    xc[0] = (xc[1] - dx[1])
    yc[0] = (yc[1] - dy[1])
    zc[0] = (zc[1] - dz[1])
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    
    ext_data = CREATE_STRUCT('x',xc, 'y',yc, 'z',zc)      

    IF (KEYWORD_SET(rho) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred density and face-centred density gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; calculate the density gradient
      dxrhof = (SHIFT(ds.rho,-1,0,0) - ds.rho) / adx
      dyrhof = (SHIFT(ds.rho,0,-1,0) - ds.rho) / ady
      dzrhof = (SHIFT(ds.rho,0,0,-1) - ds.rho) / adz

      ; enforce zero pressure gradient at maximum boundary    
      dxrhof[nx-1,*,*] = 0.0
      dyrhof[*,ny-1,*] = 0.0
      dzrhof[*,*,nz-1] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
      ;print, 'calculating cell-centred density gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxrhoc = dxrhof
      dyrhoc = dyrhof
      dzrhoc = dzrhof

      dxrhoc[1:nx-1,*,*] = (dxrhof[0:nx-2,*,*] + dxrhof[1:nx-1,*,*]) / 2.0
      dyrhoc[*,1:ny-1,*] = (dyrhof[*,0:ny-2,*] + dyrhof[*,1:ny-1,*]) / 2.0
      dzrhoc[*,*,1:nz-1] = (dzrhof[*,*,0:nz-2] + dzrhof[*,*,1:nz-1]) / 2.0
      dxrhoc[0,*,*] = 0.0
      dyrhoc[*,0,*] = 0.0
      dzrhoc[*,*,0] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('rho',ds.rho, 'dxrho',dxrhoc, 'dyrho',dyrhoc, 'dzrho',dzrhoc))
    ENDIF    
    
    IF (KEYWORD_SET(energy) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('energy',ds.energy))
    ENDIF
    
    IF (KEYWORD_SET(temperature) || KEYWORD_SET(all)) THEN BEGIN
      ;ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('temperature',ds.temperature))
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('temperature',ds.energy*(gamma-1.0d0)/2.0d0))
    ENDIF
    
    IF (KEYWORD_SET(pressure) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred pressure and face-centred pressure gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; calculate the thermal pressure
      pc = ds.rho*ds.temperature

      ; calculate the gradient of the thermal pressure
      dxpf = (SHIFT(pc,-1,0,0) - pc) / adx
      dypf = (SHIFT(pc,0,-1,0) - pc) / ady
      dzpf = (SHIFT(pc,0,0,-1) - pc) / adz

      ; enforce zero pressure gradient at maximum boundary    
      dxpf[nx-1,*,*] = 0.0
      dypf[*,ny-1,*] = 0.0
      dzpf[*,*,nz-1] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
      ;print, 'calculating cell-centred pressure gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxpc = dxpf
      dypc = dypf
      dzpc = dzpf

      dxpc[1:nx-1,*,*] = (dxpf[0:nx-2,*,*] + dxpf[1:nx-1,*,*]) / 2.0
      dypc[*,1:ny-1,*] = (dypf[*,0:ny-2,*] + dypf[*,1:ny-1,*]) / 2.0
      dzpc[*,*,1:nz-1] = (dzpf[*,*,0:nz-2] + dzpf[*,*,1:nz-1]) / 2.0
      dxpc[0,*,*] = 0.0
      dypc[*,0,*] = 0.0
      dzpc[*,*,0] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('p',pc, 'dxp',dxpc, 'dyp',dypc, 'dzp',dzpc))
    ENDIF

    IF (KEYWORD_SET(ohmic_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('ohmic_heat',ds.ohmic_heat))
    ENDIF
    
    IF (KEYWORD_SET(visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat',ds.visc_heat))
    ENDIF    
    IF (KEYWORD_SET(xy_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_xy',ds.visc_heat_xy))
    ENDIF
    IF (KEYWORD_SET(xz_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_xz',ds.visc_heat_xz))
    ENDIF
    IF (KEYWORD_SET(yz_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_yz',ds.visc_heat_yz))
    ENDIF
    IF (KEYWORD_SET(xx_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_xx',ds.visc_heat_xx))
    ENDIF
    IF (KEYWORD_SET(yy_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_yy',ds.visc_heat_yy))
    ENDIF
    IF (KEYWORD_SET(zz_visc_heat) || KEYWORD_SET(all)) THEN BEGIN
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('visc_heat_zz',ds.visc_heat_zz))
    ENDIF
    
    
    IF (KEYWORD_SET(velocity) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred velocities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      vxc = ds.vx
      vyc = ds.vy
      vzc = ds.vz

      vxc[1:nx-1,1:ny-1,1:nz-1] = (ds.vx[1:nx-1,1:ny-1,1:nz-1] + ds.vx[0:nx-2,1:ny-1,1:nz-1] + ds.vx[1:nx-1,0:ny-2,1:nz-1] + ds.vx[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.vx[0:nx-2,0:ny-2,1:nz-1] + ds.vx[0:nx-2,1:ny-1,0:nz-2] + ds.vx[1:nx-1,0:ny-2,0:nz-2] + ds.vx[0:nx-2,0:ny-2,0:nz-2]) / 8.0
      vyc[1:nx-1,1:ny-1,1:nz-1] = (ds.vy[1:nx-1,1:ny-1,1:nz-1] + ds.vy[0:nx-2,1:ny-1,1:nz-1] + ds.vy[1:nx-1,0:ny-2,1:nz-1] + ds.vy[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.vy[0:nx-2,0:ny-2,1:nz-1] + ds.vy[0:nx-2,1:ny-1,0:nz-2] + ds.vy[1:nx-1,0:ny-2,0:nz-2] + ds.vy[0:nx-2,0:ny-2,0:nz-2]) / 8.0
      vzc[1:nx-1,1:ny-1,1:nz-1] = (ds.vz[1:nx-1,1:ny-1,1:nz-1] + ds.vz[0:nx-2,1:ny-1,1:nz-1] + ds.vz[1:nx-1,0:ny-2,1:nz-1] + ds.vz[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.vz[0:nx-2,0:ny-2,1:nz-1] + ds.vz[0:nx-2,1:ny-1,0:nz-2] + ds.vz[1:nx-1,0:ny-2,0:nz-2] + ds.vz[0:nx-2,0:ny-2,0:nz-2]) / 8.0

      vxc[0,1:ny-1,1:nz-1] = (ds.vx[0,1:ny-1,1:nz-1] + ds.vx[0,0:ny-2,1:nz-1] + ds.vx[0,1:ny-1,0:nz-2] + ds.vx[0,0:ny-2,0:nz-2]) / 4.0
      vyc[0,1:ny-1,1:nz-1] = (ds.vy[0,1:ny-1,1:nz-1] + ds.vy[0,0:ny-2,1:nz-1] + ds.vy[0,1:ny-1,0:nz-2] + ds.vy[0,0:ny-2,0:nz-2]) / 4.0
      vzc[0,1:ny-1,1:nz-1] = (ds.vz[0,1:ny-1,1:nz-1] + ds.vz[0,0:ny-2,1:nz-1] + ds.vz[0,1:ny-1,0:nz-2] + ds.vz[0,0:ny-2,0:nz-2]) / 4.0
      vxc[1:nx-1,0,1:nz-1] = (ds.vx[1:nx-1,0,1:nz-1] + ds.vx[0:nx-2,0,1:nz-1] + ds.vx[1:nx-1,0,0:nz-2] + ds.vx[0:nx-2,0,0:nz-2]) / 4.0
      vyc[1:nx-1,0,1:nz-1] = (ds.vy[1:nx-1,0,1:nz-1] + ds.vy[0:nx-2,0,1:nz-1] + ds.vy[1:nx-1,0,0:nz-2] + ds.vy[0:nx-2,0,0:nz-2]) / 4.0
      vzc[1:nx-1,0,1:nz-1] = (ds.vz[1:nx-1,0,1:nz-1] + ds.vz[0:nx-2,0,1:nz-1] + ds.vz[1:nx-1,0,0:nz-2] + ds.vz[0:nx-2,0,0:nz-2]) / 4.0
      vxc[1:nx-1,1:ny-1,0] = (ds.vx[1:nx-1,1:ny-1,0] + ds.vx[0:nx-2,1:ny-1,0] + ds.vx[1:nx-1,0:ny-2,0] + ds.vx[0:nx-2,0:ny-2,0]) / 4.0
      vyc[1:nx-1,1:ny-1,0] = (ds.vy[1:nx-1,1:ny-1,0] + ds.vy[0:nx-2,1:ny-1,0] + ds.vy[1:nx-1,0:ny-2,0] + ds.vy[0:nx-2,0:ny-2,0]) / 4.0
      vzc[1:nx-1,1:ny-1,0] = (ds.vz[1:nx-1,1:ny-1,0] + ds.vz[0:nx-2,1:ny-1,0] + ds.vz[1:nx-1,0:ny-2,0] + ds.vz[0:nx-2,0:ny-2,0]) / 4.0
    
      vxc[0,0,1:nz-1] = (ds.vx[0,0,1:nz-1] + ds.vx[0,0,0:nz-2]) / 2.0
      vyc[0,0,1:nz-1] = (ds.vy[0,0,1:nz-1] + ds.vy[0,0,0:nz-2]) / 2.0
      vzc[0,0,1:nz-1] = (ds.vz[0,0,1:nz-1] + ds.vz[0,0,0:nz-2]) / 2.0
      vxc[0,1:ny-1,0] = (ds.vx[0,1:ny-1,0] + ds.vx[0,0:ny-2,0]) / 2.0
      vyc[0,1:ny-1,0] = (ds.vy[0,1:ny-1,0] + ds.vy[0,0:ny-2,0]) / 2.0
      vzc[0,1:ny-1,0] = (ds.vz[0,1:ny-1,0] + ds.vz[0,0:ny-2,0]) / 2.0
      vxc[1:nx-1,0,0] = (ds.vx[1:nx-1,0,0] + ds.vx[0:nx-2,0,0]) / 2.0
      vyc[1:nx-1,0,0] = (ds.vy[1:nx-1,0,0] + ds.vy[0:nx-2,0,0]) / 2.0
      vzc[1:nx-1,0,0] = (ds.vz[1:nx-1,0,0] + ds.vz[0:nx-2,0,0]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('vx',vxc,'vy',vyc,'vz',vzc))
    ENDIF
    
    IF (KEYWORD_SET(vorticity) || KEYWORD_SET(all)) THEN BEGIN    
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
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ;print, 'calculating edge-centred vorticities (using face-centred velocity components)...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      wxe = ((SHIFT(vzf,0,-1,0) - vzf) / ady) - ((SHIFT(vyf,0,0,-1) - vyf) / adz)
      wye = ((SHIFT(vxf,0,0,-1) - vxf) / adz) - ((SHIFT(vzf,-1,0,0) - vzf) / adx)
      wze = ((SHIFT(vyf,-1,0,0) - vyf) / adx) - ((SHIFT(vxf,0,-1,0) - vxf) / ady)

      ; enforce zero vorticity since velocities are zero at boundary
      wxe[nx-1,*,*] = zero
      wye[*,ny-1,*] = zero
      wze[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating cell-centred vorticities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      wxc = wxe
      wyc = wye
      wzc = wze

      wxc[*,1:ny-1,1:nz-1] = (wxe[*,1:ny-1,1:nz-1] + wxe[*,0:ny-2,1:nz-1] + wxe[*,1:ny-1,0:nz-2] + wxe[*,0:ny-2,0:nz-2]) / 4.0
      wyc[1:nx-1,*,1:nz-1] = (wye[1:nx-1,*,1:nz-1] + wye[0:nx-2,*,1:nz-1] + wye[1:nx-1,*,0:nz-2] + wye[0:nx-2,*,0:nz-2]) / 4.0
      wzc[1:nx-1,1:ny-1,*] = (wze[1:nx-1,1:ny-1,*] + wze[0:nx-2,1:ny-1,*] + wze[1:nx-1,0:ny-2,*] + wze[0:nx-2,0:ny-2,*]) / 4.0

      wxc[*,0,1:nz-1] = (wxe[*,0,1:nz-1] + wxe[*,0,0:nz-2]) / 2.0
      wyc[0,*,1:nz-1] = (wye[0,*,1:nz-1] + wye[0,*,0:nz-2]) / 2.0
      wzc[0,1:ny-1,*] = (wze[0,1:ny-1,*] + wze[0,0:ny-2,*]) / 2.0

      wxc[*,1:ny-1,0] = (wxe[*,1:ny-1,0] + wxe[*,0:ny-2,0]) / 2.0 
      wyc[1:nx-1,*,0] = (wye[1:nx-1,*,0] + wye[0:nx-2,*,0]) / 2.0   
      wzc[1:nx-1,0,*] = (wze[1:nx-1,0,*] + wze[0:nx-2,0,*]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('wx',wxc,'wy',wyc,'wz',wzc))
    ENDIF

    IF (KEYWORD_SET(bfield) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred magnetic field components...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      bxc = (SHIFT(ds.bx,-1,0,0) + ds.bx) / 2.0
      byc = (SHIFT(ds.by,0,-1,0) + ds.by) / 2.0
      bzc = (SHIFT(ds.bz,0,0,-1) + ds.bz) / 2.0

      ; enforce zero magnetic gradient at minimum boundary
      bxc[nx-1,*,*] = bxc[nx-2,*,*]
      byc[*,ny-1,*] = byc[*,ny-2,*]
      bzc[*,*,nz-1] = bzc[*,*,nz-2]
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('bx',bxc,'by',byc,'bz',bzc))
    ENDIF

    IF (KEYWORD_SET(current) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating edge-centred currents...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      jxe = ((SHIFT(ds.bz,0,-1,0) - ds.bz) / ady) - ((SHIFT(ds.by,0,0,-1) - ds.by) / adz)
      jye = ((SHIFT(ds.bx,0,0,-1) - ds.bx) / adz) - ((SHIFT(ds.bz,-1,0,0) - ds.bz) / adx)
      jze = ((SHIFT(ds.by,-1,0,0) - ds.by) / adx) - ((SHIFT(ds.bx,0,-1,0) - ds.bx) / ady)

      ; enforce zero current since magnetic field gradient is zero at minimum boundary
      jxe[nx-1,*,*] = zero
      jye[*,ny-1,*] = zero
      jze[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;        
    
      ;print, 'calculating cell-centred currents...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      jxc = jxe
      jyc = jye
      jzc = jze

      jxc[*,1:ny-1,1:nz-1] = (jxe[*,1:ny-1,1:nz-1] + jxe[*,0:ny-2,1:nz-1] + jxe[*,1:ny-1,0:nz-2] + jxe[*,0:ny-2,0:nz-2]) / 4.0
      jyc[1:nx-1,*,1:nz-1] = (jye[1:nx-1,*,1:nz-1] + jye[0:nx-2,*,1:nz-1] + jye[1:nx-1,*,0:nz-2] + jye[0:nx-2,*,0:nz-2]) / 4.0
      jzc[1:nx-1,1:ny-1,*] = (jze[1:nx-1,1:ny-1,*] + jze[0:nx-2,1:ny-1,*] + jze[1:nx-1,0:ny-2,*] + jze[0:nx-2,0:ny-2,*]) / 4.0

      jxc[*,0,1:nz-1] = (jxe[*,0,1:nz-1] + jxe[*,0,0:nz-2]) / 2.0
      jyc[0,*,1:nz-1] = (jye[0,*,1:nz-1] + jye[0,*,0:nz-2]) / 2.0
      jzc[0,1:ny-1,*] = (jze[0,1:ny-1,*] + jze[0,0:ny-2,*]) / 2.0

      jxc[*,1:ny-1,0] = (jxe[*,1:ny-1,0] + jxe[*,0:ny-2,0]) / 2.0 
      jyc[1:nx-1,*,0] = (jye[1:nx-1,*,0] + jye[0:nx-2,*,0]) / 2.0   
      jzc[1:nx-1,0,*] = (jze[1:nx-1,0,*] + jze[0:nx-2,0,*]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('jx',jxc,'jy',jyc,'jz',jzc))
    ENDIF
    
    IF ((KEYWORD_SET(current) && KEYWORD_SET(bfield) && KEYWORD_SET(lorentz)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating cell-centred Lorentz forces...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      fxc = jyc*bzc - jzc*byc
      fyc = jzc*bxc - jxc*bzc
      fzc = jxc*byc - jyc*bxc
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('fx',fxc, 'fy',fyc, 'fz',fzc))
    ENDIF
    
    
    
    IF (KEYWORD_SET(ffield) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred raw lorentz forces...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      fxc = ds.fx
      fyc = ds.fy
      fzc = ds.fz

      fxc[1:nx-1,1:ny-1,1:nz-1] = (ds.fx[1:nx-1,1:ny-1,1:nz-1] + ds.fx[0:nx-2,1:ny-1,1:nz-1] + ds.fx[1:nx-1,0:ny-2,1:nz-1] + ds.fx[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.fx[0:nx-2,0:ny-2,1:nz-1] + ds.fx[0:nx-2,1:ny-1,0:nz-2] + ds.fx[1:nx-1,0:ny-2,0:nz-2] + ds.fx[0:nx-2,0:ny-2,0:nz-2]) / 8.0
      fyc[1:nx-1,1:ny-1,1:nz-1] = (ds.fy[1:nx-1,1:ny-1,1:nz-1] + ds.fy[0:nx-2,1:ny-1,1:nz-1] + ds.fy[1:nx-1,0:ny-2,1:nz-1] + ds.fy[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.fy[0:nx-2,0:ny-2,1:nz-1] + ds.fy[0:nx-2,1:ny-1,0:nz-2] + ds.fy[1:nx-1,0:ny-2,0:nz-2] + ds.fy[0:nx-2,0:ny-2,0:nz-2]) / 8.0
      fzc[1:nx-1,1:ny-1,1:nz-1] = (ds.fz[1:nx-1,1:ny-1,1:nz-1] + ds.fz[0:nx-2,1:ny-1,1:nz-1] + ds.fz[1:nx-1,0:ny-2,1:nz-1] + ds.fz[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.fz[0:nx-2,0:ny-2,1:nz-1] + ds.fz[0:nx-2,1:ny-1,0:nz-2] + ds.fz[1:nx-1,0:ny-2,0:nz-2] + ds.fz[0:nx-2,0:ny-2,0:nz-2]) / 8.0

      fxc[0,1:ny-1,1:nz-1] = (ds.fx[0,1:ny-1,1:nz-1] + ds.fx[0,0:ny-2,1:nz-1] + ds.fx[0,1:ny-1,0:nz-2] + ds.fx[0,0:ny-2,0:nz-2]) / 4.0
      fyc[0,1:ny-1,1:nz-1] = (ds.fy[0,1:ny-1,1:nz-1] + ds.fy[0,0:ny-2,1:nz-1] + ds.fy[0,1:ny-1,0:nz-2] + ds.fy[0,0:ny-2,0:nz-2]) / 4.0
      fzc[0,1:ny-1,1:nz-1] = (ds.fz[0,1:ny-1,1:nz-1] + ds.fz[0,0:ny-2,1:nz-1] + ds.fz[0,1:ny-1,0:nz-2] + ds.fz[0,0:ny-2,0:nz-2]) / 4.0
      fxc[1:nx-1,0,1:nz-1] = (ds.fx[1:nx-1,0,1:nz-1] + ds.fx[0:nx-2,0,1:nz-1] + ds.fx[1:nx-1,0,0:nz-2] + ds.fx[0:nx-2,0,0:nz-2]) / 4.0
      fyc[1:nx-1,0,1:nz-1] = (ds.fy[1:nx-1,0,1:nz-1] + ds.fy[0:nx-2,0,1:nz-1] + ds.fy[1:nx-1,0,0:nz-2] + ds.fy[0:nx-2,0,0:nz-2]) / 4.0
      fzc[1:nx-1,0,1:nz-1] = (ds.fz[1:nx-1,0,1:nz-1] + ds.fz[0:nx-2,0,1:nz-1] + ds.fz[1:nx-1,0,0:nz-2] + ds.fz[0:nx-2,0,0:nz-2]) / 4.0
      fxc[1:nx-1,1:ny-1,0] = (ds.fx[1:nx-1,1:ny-1,0] + ds.fx[0:nx-2,1:ny-1,0] + ds.fx[1:nx-1,0:ny-2,0] + ds.fx[0:nx-2,0:ny-2,0]) / 4.0
      fyc[1:nx-1,1:ny-1,0] = (ds.fy[1:nx-1,1:ny-1,0] + ds.fy[0:nx-2,1:ny-1,0] + ds.fy[1:nx-1,0:ny-2,0] + ds.fy[0:nx-2,0:ny-2,0]) / 4.0
      fzc[1:nx-1,1:ny-1,0] = (ds.fz[1:nx-1,1:ny-1,0] + ds.fz[0:nx-2,1:ny-1,0] + ds.fz[1:nx-1,0:ny-2,0] + ds.fz[0:nx-2,0:ny-2,0]) / 4.0
    
      fxc[0,0,1:nz-1] = (ds.fx[0,0,1:nz-1] + ds.fx[0,0,0:nz-2]) / 2.0
      fyc[0,0,1:nz-1] = (ds.fy[0,0,1:nz-1] + ds.fy[0,0,0:nz-2]) / 2.0
      fzc[0,0,1:nz-1] = (ds.fz[0,0,1:nz-1] + ds.fz[0,0,0:nz-2]) / 2.0
      fxc[0,1:ny-1,0] = (ds.fx[0,1:ny-1,0] + ds.fx[0,0:ny-2,0]) / 2.0
      fyc[0,1:ny-1,0] = (ds.fy[0,1:ny-1,0] + ds.fy[0,0:ny-2,0]) / 2.0
      fzc[0,1:ny-1,0] = (ds.fz[0,1:ny-1,0] + ds.fz[0,0:ny-2,0]) / 2.0
      fxc[1:nx-1,0,0] = (ds.fx[1:nx-1,0,0] + ds.fx[0:nx-2,0,0]) / 2.0
      fyc[1:nx-1,0,0] = (ds.fy[1:nx-1,0,0] + ds.fy[0:nx-2,0,0]) / 2.0
      fzc[1:nx-1,0,0] = (ds.fz[1:nx-1,0,0] + ds.fz[0:nx-2,0,0]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('fx',fxc,'fy',fyc,'fz',fzc))
    ENDIF
    
    

    IF (KEYWORD_SET(resistivity) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred resistivities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;             
      etac = ds.eta
      
      etac[1:nx-1,1:ny-1,1:nz-1] = (ds.eta[1:nx-1,1:ny-1,1:nz-1] + ds.eta[0:nx-2,1:ny-1,1:nz-1] + ds.eta[1:nx-1,0:ny-2,1:nz-1] + ds.eta[1:nx-1,1:ny-1,0:nz-2]$
                                 + ds.eta[0:nx-2,0:ny-2,1:nz-1] + ds.eta[0:nx-2,1:ny-1,0:nz-2] + ds.eta[1:nx-1,0:ny-2,0:nz-2] + ds.eta[0:nx-2,0:ny-2,0:nz-2]) / 8.0
      
      etac[0,1:ny-1,1:nz-1] = (ds.eta[0,1:ny-1,1:nz-1] + ds.eta[0,0:ny-2,1:nz-1] + ds.eta[0,1:ny-1,0:nz-2] + ds.eta[0,0:ny-2,0:nz-2]) / 4.0
      etac[1:nx-1,0,1:nz-1] = (ds.eta[1:nx-1,0,1:nz-1] + ds.eta[0:nx-2,0,1:nz-1] + ds.eta[1:nx-1,0,0:nz-2] + ds.eta[0:nx-2,0,0:nz-2]) / 4.0
      etac[1:nx-1,1:ny-1,0] = (ds.eta[1:nx-1,1:ny-1,0] + ds.eta[0:nx-2,1:ny-1,0] + ds.eta[1:nx-1,0:ny-2,0] + ds.eta[0:nx-2,0:ny-2,0]) / 4.0
      
      etac[0,0,1:nz-1] = (ds.eta[0,0,1:nz-1] + ds.eta[0,0,0:nz-2]) / 2.0
      etac[0,1:ny-1,0] = (ds.eta[0,1:ny-1,0] + ds.eta[0,0:ny-2,0]) / 2.0
      etac[1:nx-1,0,0] = (ds.eta[1:nx-1,0,0] + ds.eta[0:nx-2,0,0]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('eta',etac))
    ENDIF
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    IF ((KEYWORD_SET(resistivity) && KEYWORD_SET(current) && KEYWORD_SET(bfield) && KEYWORD_SET(efield)) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating cell-centred electric field components (mu0=1)...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      exc = etac*jxc - (vyc*bzc - vzc*byc)
      eyc = etac*jyc - (vzc*bxc - vxc*bzc)
      ezc = etac*jzc - (vxc*byc - vyc*bxc)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('ex',exc,'ey',eyc,'ez',ezc))
    ENDIF


    IF ((KEYWORD_SET(current) && KEYWORD_SET(bfield) && KEYWORD_SET(alpha)) || KEYWORD_SET(all)) THEN BEGIN      
      ;print, 'calculating cell-centred alpha...'
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      alphac = jxc
      alphac[*,*,*] = 0.0
      bmag_nz_i = WHERE((bxc^2 + byc^2 + bzc^2) GT 0.0, cnt)
      IF cnt GT 0 THEN BEGIN
        x_i = bmag_nz_i mod nx
        y_i = (bmag_nz_i/nx) mod ny
        z_i = bmag_nz_i/(nx*ny)
        alphac[x_i,y_i,z_i] = (bxc[x_i,y_i,z_i]*jxc[x_i,y_i,z_i] + byc[x_i,y_i,z_i]*jyc[x_i,y_i,z_i] + bzc[x_i,y_i,z_i]*jzc[x_i,y_i,z_i])$
                             /(bxc[x_i,y_i,z_i]^2 + byc[x_i,y_i,z_i]^2 + bzc[x_i,y_i,z_i]^2)
      ENDIF
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ext_data = CREATE_STRUCT(ext_data, 'alpha',alphac, 'nff',nffc, 'dtbmag',dtbmagc)        
    ENDIF

    
    IF ((KEYWORD_SET(current) && KEYWORD_SET(resistivity) && KEYWORD_SET(bfield) && KEYWORD_SET(efield) && KEYWORD_SET(dtW)) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'calculating face-centred electric field components (mu0=1)...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      exf = exc
      eyf = eyc
      ezf = ezc

      exf[0:nx-2,*,*] = (exc[0:nx-2,*,*] + exc[1:nx-1,*,*]) / 2.0
      eyf[*,0:ny-2,*] = (eyc[*,0:ny-2,*] + eyc[*,1:ny-1,*]) / 2.0
      ezf[*,*,0:nz-2] = (ezc[*,*,0:nz-2] + ezc[*,*,1:nz-1]) / 2.0

      exf[nx-1,*,*] = 0.0
      eyf[*,ny-1,*] = 0.0
      ezf[*,*,nz-1] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating edge-centred curl of electric field...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      curlexe = ((SHIFT(ezf,0,-1,0) - ezf) / ady) - ((SHIFT(eyf,0,0,-1) - eyf) / adz)
      curleye = ((SHIFT(exf,0,0,-1) - exf) / adz) - ((SHIFT(ezf,-1,0,0) - ezf) / adx)
      curleze = ((SHIFT(eyf,-1,0,0) - eyf) / adx) - ((SHIFT(exf,0,-1,0) - exf) / ady)

      ; enforce zero values since electric field gradient is zero at boundaries
      curlexe[nx-1,*,*] = zero
      curleye[*,ny-1,*] = zero
      curleze[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

      ;print, 'calculating cell-centred curl of electric field...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      curlexc = curlexe
      curleyc = curleye
      curlezc = curleze

      curlexc[*,1:ny-1,1:nz-1] = (curlexe[*,1:ny-1,1:nz-1] + curlexe[*,0:ny-2,1:nz-1] + curlexe[*,1:ny-1,0:nz-2] + curlexe[*,0:ny-2,0:nz-2]) / 4.0
      curleyc[1:nx-1,*,1:nz-1] = (curleye[1:nx-1,*,1:nz-1] + curleye[0:nx-2,*,1:nz-1] + curleye[1:nx-1,*,0:nz-2] + curleye[0:nx-2,*,0:nz-2]) / 4.0
      curlezc[1:nx-1,1:ny-1,*] = (curleze[1:nx-1,1:ny-1,*] + curleze[0:nx-2,1:ny-1,*] + curleze[1:nx-1,0:ny-2,*] + curleze[0:nx-2,0:ny-2,*]) / 4.0

      curlexc[*,0,1:nz-1] = (curlexe[*,0,1:nz-1] + curlexe[*,0,0:nz-2]) / 2.0
      curleyc[0,*,1:nz-1] = (curleye[0,*,1:nz-1] + curleye[0,*,0:nz-2]) / 2.0
      curlezc[0,1:ny-1,*] = (curleze[0,1:ny-1,*] + curleze[0,0:ny-2,*]) / 2.0

      curlexc[*,1:ny-1,0] = (curlexe[*,1:ny-1,0] + curlexe[*,0:ny-2,0]) / 2.0 
      curleyc[1:nx-1,*,0] = (curleye[1:nx-1,*,0] + curleye[0:nx-2,*,0]) / 2.0   
      curlezc[1:nx-1,0,*] = (curleze[1:nx-1,0,*] + curleze[0:nx-2,0,*]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating cell-centred time derivative of the magnetic field strength...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dtWxc = -(bxc*curlexc)
      dtWyc = -(byc*curleyc)
      dtWzc = -(bzc*curlezc)
      ;dtWc = -(bxc*curlexc + byc*curleyc + bzc*curlezc)

      ; enforce zero values since electric field gradient is zero at boundaries
      dtWxc[nx-1,*,*] = zero
      dtWxc[*,ny-1,*] = zero
      dtWxc[*,*,nz-1] = zero
      dtWyc[nx-1,*,*] = zero
      dtWyc[*,ny-1,*] = zero
      dtWyc[*,*,nz-1] = zero
      dtWzc[nx-1,*,*] = zero
      dtWzc[*,ny-1,*] = zero
      dtWzc[*,*,nz-1] = zero
      ;dtWc[nx-1,*,*] = zero
      ;dtWc[*,ny-1,*] = zero
      ;dtWc[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       

      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('dtWx',dtWxc,'dtWy',dtWyc,'dtWz',dtWzc))
      ;ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('dtW',dtWc)) 
    ENDIF 
    
    
    
                                         

  ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data structure>"
    ext_data = 0
  ENDELSE  

  RETURN, ext_data

END
