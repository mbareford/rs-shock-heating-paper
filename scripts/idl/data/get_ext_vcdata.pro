; Return an extended data set for a particular snapshot.
;
; The data set comprises the following.
;
; Cell boundary coordinates.
; Vertex-centred magnetic field components.
; Vertex-centred velocity components.
; Vertex-centred currents, resistivities and vorticities (calculated from face-centred velocities).
; Vertex-centred Lorentz forces.
;
; Vertex-centred electric field components calculated from vertex-centred currents and resistivities,
; Vertex-centred magnetic field components, and vertex-centred velocity components.
;
; Alpha values calculated from vertex-centred currents and vertex-centred magnetic field components.
;
; ds: the basic data set for some snapshot
; jcrit: the current threshold above which anomalous resistivity is applied
; eta0: the background resistivity
; eta1: the anomalous resistivity
FUNCTION get_ext_vcdata, ds, jcrit, eta0, eta1, density=density, pressure=pressure, $
                                                velocity=velocity, vorticity=vorticity, $
                                                bfield=bfield, efield=efield, resistivity=resistivity, current=current, $
                                                lorentz=lorentz, alpha=alpha, dtW=dtW, all=all

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
    ext_data = CREATE_STRUCT('x',ds.grid.x, 'y',ds.grid.y, 'z',ds.grid.z)    


    IF (KEYWORD_SET(density) || KEYWORD_SET(all)) THEN BEGIN
      
      ;print, 'calculating the vertex-centred density...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      rhov = ds.rho
      
      rhov[0:nx-2,0:ny-2,0:nz-2] = (  ds.rho[0:nx-2,0:ny-2,0:nz-2] + ds.rho[1:nx-1,0:ny-2,0:nz-2] + ds.rho[0:nx-2,1:ny-1,0:nz-2] + ds.rho[1:nx-1,1:ny-1,0:nz-2] $
                                    + ds.rho[0:nx-2,0:ny-2,1:nz-1] + ds.rho[1:nx-1,0:ny-2,1:nz-1] + ds.rho[0:nx-2,1:ny-1,1:nz-1] + ds.rho[1:nx-1,1:ny-1,1:nz-1]) / 8.0 

      ; enforce zero gradient at maximum boundaries
      rhov[nx-1,0:ny-2,0:nz-2] = (ds.rho[nx-1,0:ny-2,0:nz-2] + ds.rho[nx-1,1:ny-1,0:nz-2] + ds.rho[nx-1,0:ny-2,1:nz-1] + ds.rho[nx-1,1:ny-1,1:nz-1]) / 4.0
      rhov[0:nx-2,ny-1,0:nz-2] = (ds.rho[0:nx-2,ny-1,0:nz-2] + ds.rho[1:nx-1,ny-1,0:nz-2] + ds.rho[0:nx-2,ny-1,1:nz-1] + ds.rho[1:nx-1,ny-1,1:nz-1]) / 4.0
      rhov[0:nx-2,0:ny-2,nz-1] = (ds.rho[0:nx-2,0:ny-2,nz-1] + ds.rho[1:nx-1,0:ny-2,nz-1] + ds.rho[0:nx-2,1:ny-1,nz-1] + ds.rho[1:nx-1,1:ny-1,nz-1]) / 4.0      
      rhov[nx-1,ny-1,0:nz-2]   = (ds.rho[nx-1,ny-1,0:nz-2] + ds.rho[nx-1,ny-1,1:nz-1]) / 2.0
      rhov[nx-1,0:ny-2,nz-1]   = (ds.rho[nx-1,0:ny-2,nz-1] + ds.rho[nx-1,1:ny-1,nz-1]) / 2.0
      rhov[0:nx-2,ny-1,nz-1]   = (ds.rho[0:nx-2,ny-1,nz-1] + ds.rho[1:nx-1,ny-1,nz-1]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating face-centred density gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxrhof = (SHIFT(ds.rho,-1,0,0) - ds.rho) / adx
      dyrhof = (SHIFT(ds.rho,0,-1,0) - ds.rho) / ady
      dzrhof = (SHIFT(ds.rho,0,0,-1) - ds.rho) / adz

      ; enforce zero gradient at maximum boundaries
      dxrhof[nx-1,*,*] = 0.0
      dyrhof[*,ny-1,*] = 0.0
      dzrhof[*,*,nz-1] = 0.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating vertex-centred density gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxrhov = dxrhof
      dyrhov = dyrhof
      dzrhov = dzrhof

      dxrhov[0:nx-1,0:ny-2,0:nz-2] = (dxrhof[0:nx-1,0:ny-2,0:nz-2] + dxrhof[0:nx-1,1:ny-1,0:nz-2] + dxrhof[0:nx-1,0:ny-2,1:nz-1] + dxrhof[0:nx-1,1:ny-1,1:nz-1]) / 4.0
      dyrhov[0:nx-2,0:ny-1,0:nz-2] = (dyrhof[0:nx-2,0:ny-1,0:nz-2] + dyrhof[1:nx-1,0:ny-1,0:nz-2] + dyrhof[0:nx-2,0:ny-1,1:nz-1] + dyrhof[1:nx-1,0:ny-1,1:nz-1]) / 4.0
      dzrhov[0:nx-2,0:ny-2,0:nz-1] = (dzrhof[0:nx-2,0:ny-2,0:nz-1] + dzrhof[1:nx-1,0:ny-2,0:nz-1] + dzrhof[0:nx-2,1:ny-1,0:nz-1] + dzrhof[1:nx-1,1:ny-1,0:nz-1]) / 4.0

      dxrhov[0:nx-1,ny-1,0:nz-2] = (dxrhof[0:nx-1,ny-1,0:nz-2] + dxrhof[0:nx-1,ny-1,1:nz-1]) / 2.0
      dxrhov[0:nx-1,0:ny-2,nz-1] = (dxrhof[0:nx-1,0:ny-2,nz-1] + dxrhof[0:nx-1,1:ny-1,nz-1]) / 2.0  
      dyrhov[0:nx-2,0:ny-1,nz-1] = (dyrhof[0:nx-2,0:ny-1,nz-1] + dyrhof[1:nx-1,0:ny-1,nz-1]) / 2.0      
      dyrhov[nx-1,0:ny-1,0:nz-2] = (dyrhof[nx-1,0:ny-1,0:nz-2] + dyrhof[nx-1,0:ny-1,1:nz-1]) / 2.0
      dzrhov[0:nx-2,ny-1,0:nz-1] = (dzrhof[0:nx-2,ny-1,0:nz-1] + dzrhof[1:nx-1,ny-1,0:nz-1]) / 2.0  
      dzrhov[nx-1,0:ny-2,0:nz-1] = (dzrhof[nx-1,0:ny-2,0:nz-1] + dzrhof[nx-1,1:ny-1,0:nz-1]) / 2.0     
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('rho',rhov, 'dxrho',dxrhov, 'dyrho',dyrhov, 'dzrho',dzrhov))
    ENDIF


    IF (KEYWORD_SET(pressure) || KEYWORD_SET(all)) THEN BEGIN

      ;print, 'calculating the vertex-centred pressure...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      pc = ds.rho*ds.temperature
      pv = pc
      
      pv[0:nx-2,0:ny-2,0:nz-2] = (  pc[0:nx-2,0:ny-2,0:nz-2] + pc[1:nx-1,0:ny-2,0:nz-2] + pc[0:nx-2,1:ny-1,0:nz-2] + pc[1:nx-1,1:ny-1,0:nz-2] $
                                    + pc[0:nx-2,0:ny-2,1:nz-1] + pc[1:nx-1,0:ny-2,1:nz-1] + pc[0:nx-2,1:ny-1,1:nz-1] + pc[1:nx-1,1:ny-1,1:nz-1]) / 8.0 

      ; enforce zero gradient at maximum boundaries
      pv[nx-1,0:ny-2,0:nz-2] = (pc[nx-1,0:ny-2,0:nz-2] + pc[nx-1,1:ny-1,0:nz-2] + pc[nx-1,0:ny-2,1:nz-1] + pc[nx-1,1:ny-1,1:nz-1]) / 4.0
      pv[0:nx-2,ny-1,0:nz-2] = (pc[0:nx-2,ny-1,0:nz-2] + pc[1:nx-1,ny-1,0:nz-2] + pc[0:nx-2,ny-1,1:nz-1] + pc[1:nx-1,ny-1,1:nz-1]) / 4.0
      pv[0:nx-2,0:ny-2,nz-1] = (pc[0:nx-2,0:ny-2,nz-1] + pc[1:nx-1,0:ny-2,nz-1] + pc[0:nx-2,1:ny-1,nz-1] + pc[1:nx-1,1:ny-1,nz-1]) / 4.0      
      pv[nx-1,ny-1,0:nz-2]   = (pc[nx-1,ny-1,0:nz-2] + pc[nx-1,ny-1,1:nz-1]) / 2.0
      pv[nx-1,0:ny-2,nz-1]   = (pc[nx-1,0:ny-2,nz-1] + pc[nx-1,1:ny-1,nz-1]) / 2.0
      pv[0:nx-2,ny-1,nz-1]   = (pc[0:nx-2,ny-1,nz-1] + pc[1:nx-1,ny-1,nz-1]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ;print, 'calculating the face-centred pressure gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxpf = (SHIFT(pc,-1,0,0) - pc) / adx
      dypf = (SHIFT(pc,0,-1,0) - pc) / ady
      dzpf = (SHIFT(pc,0,0,-1) - pc) / adz

      ; enforce zero gradient at maximum boundaries
      dxpf[nx-1,*,*] = 0.0
      dypf[*,ny-1,*] = 0.0
      dzpf[*,*,nz-1] = 0.0     
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;print, 'calculating vertex-centred pressure gradients...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dxpv = dxpf
      dypv = dypf
      dzpv = dzpf

      dxpv[0:nx-1,0:ny-2,0:nz-2] = (dxpf[0:nx-1,0:ny-2,0:nz-2] + dxpf[0:nx-1,1:ny-1,0:nz-2] + dxpf[0:nx-1,0:ny-2,1:nz-1] + dxpf[0:nx-1,1:ny-1,1:nz-1]) / 4.0
      dypv[0:nx-2,0:ny-1,0:nz-2] = (dypf[0:nx-2,0:ny-1,0:nz-2] + dypf[1:nx-1,0:ny-1,0:nz-2] + dypf[0:nx-2,0:ny-1,1:nz-1] + dypf[1:nx-1,0:ny-1,1:nz-1]) / 4.0
      dzpv[0:nx-2,0:ny-2,0:nz-1] = (dzpf[0:nx-2,0:ny-2,0:nz-1] + dzpf[1:nx-1,0:ny-2,0:nz-1] + dzpf[0:nx-2,1:ny-1,0:nz-1] + dzpf[1:nx-1,1:ny-1,0:nz-1]) / 4.0

      dxpv[0:nx-1,ny-1,0:nz-2] = (dxpf[0:nx-1,ny-1,0:nz-2] + dxpf[0:nx-1,ny-1,1:nz-1]) / 2.0
      dxpv[0:nx-1,0:ny-2,nz-1] = (dxpf[0:nx-1,0:ny-2,nz-1] + dxpf[0:nx-1,1:ny-1,nz-1]) / 2.0  
      dypv[0:nx-2,0:ny-1,nz-1] = (dypf[0:nx-2,0:ny-1,nz-1] + dypf[1:nx-1,0:ny-1,nz-1]) / 2.0      
      dypv[nx-1,0:ny-1,0:nz-2] = (dypf[nx-1,0:ny-1,0:nz-2] + dypf[nx-1,0:ny-1,1:nz-1]) / 2.0
      dzpv[0:nx-2,ny-1,0:nz-1] = (dzpf[0:nx-2,ny-1,0:nz-1] + dzpf[1:nx-1,ny-1,0:nz-1]) / 2.0  
      dzpv[nx-1,0:ny-2,0:nz-1] = (dzpf[nx-1,0:ny-2,0:nz-1] + dzpf[nx-1,1:ny-1,0:nz-1]) / 2.0 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('p',pv, 'dxp',dxpv, 'dyp',dypv, 'dzp',dzpv))
    ENDIF


    IF (KEYWORD_SET(velocity) || KEYWORD_SET(all)) THEN BEGIN
      ;print, 'using vertex-centred velocities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('vx',ds.vx, 'vy',ds.vy, 'vz',ds.vz))      
    ENDIF


    IF (KEYWORD_SET(bfield) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating vertex-centred magnetic field components...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      bxv = ds.bx
      byv = ds.by
      bzv = ds.bz

      bxv[0:nx-1,0:ny-2,0:nz-2] = (ds.bx[0:nx-1,0:ny-2,0:nz-2] + ds.bx[0:nx-1,1:ny-1,0:nz-2] + ds.bx[0:nx-1,0:ny-2,1:nz-1] + ds.bx[0:nx-1,1:ny-1,1:nz-1]) / 4.0
      byv[0:nx-2,0:ny-1,0:nz-2] = (ds.by[0:nx-2,0:ny-1,0:nz-2] + ds.by[1:nx-1,0:ny-1,0:nz-2] + ds.by[0:nx-2,0:ny-1,1:nz-1] + ds.by[1:nx-1,0:ny-1,1:nz-1]) / 4.0
      bzv[0:nx-2,0:ny-2,0:nz-1] = (ds.bz[0:nx-2,0:ny-2,0:nz-1] + ds.bz[1:nx-1,0:ny-2,0:nz-1] + ds.bz[0:nx-2,1:ny-1,0:nz-1] + ds.bz[1:nx-1,1:ny-1,0:nz-1]) / 4.0

      bxv[0:nx-1,ny-1,0:nz-2] = (ds.bx[0:nx-1,ny-1,0:nz-2] + ds.bx[0:nx-1,ny-1,1:nz-1]) / 2.0
      bxv[0:nx-1,0:ny-2,nz-1] = (ds.bx[0:nx-1,0:ny-2,nz-1] + ds.bx[0:nx-1,1:ny-1,nz-1]) / 2.0  
      byv[0:nx-2,0:ny-1,nz-1] = (ds.by[0:nx-2,0:ny-1,nz-1] + ds.by[1:nx-1,0:ny-1,nz-1]) / 2.0      
      byv[nx-1,0:ny-1,0:nz-2] = (ds.by[nx-1,0:ny-1,0:nz-2] + ds.by[nx-1,0:ny-1,1:nz-1]) / 2.0
      bzv[0:nx-2,ny-1,0:nz-1] = (ds.bz[0:nx-2,ny-1,0:nz-1] + ds.bz[1:nx-1,ny-1,0:nz-1]) / 2.0  
      bzv[nx-1,0:ny-2,0:nz-1] = (ds.bz[nx-1,0:ny-2,0:nz-1] + ds.bz[nx-1,1:ny-1,0:nz-1]) / 2.0
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('bx',bxv, 'by',byv, 'bz',bzv))
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

      ;print, 'calculating vertex-centred vorticities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      wxv = wxe
      wyv = wye
      wzv = wze

      wxv[0:nx-2,0:ny-1,0:nz-1] = (wxe[0:nx-2,0:ny-1,0:nz-1] + wxe[1:nx-1,0:ny-1,0:nz-1]) / 2.0
      wyv[0:nx-1,0:ny-2,0:nz-1] = (wye[0:nx-1,0:ny-2,0:nz-1] + wye[0:nx-1,1:ny-1,0:nz-1]) / 2.0
      wzv[0:nx-1,0:ny-1,0:nz-2] = (wze[0:nx-1,0:ny-1,0:nz-2] + wze[0:nx-1,0:ny-1,1:nz-1]) / 2.0    
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('wx',wxv, 'wy',wyv, 'wz',wzv))
    ENDIF


    IF (KEYWORD_SET(current) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating edge-centred currents (using face-centred magnetic field components)...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      jxe = ((SHIFT(ds.bz,0,-1,0) - ds.bz) / ady) - ((SHIFT(ds.by,0,0,-1) - ds.by) / adz)
      jye = ((SHIFT(ds.bx,0,0,-1) - ds.bx) / adz) - ((SHIFT(ds.bz,-1,0,0) - ds.bz) / adx)
      jze = ((SHIFT(ds.by,-1,0,0) - ds.by) / adx) - ((SHIFT(ds.bx,0,-1,0) - ds.bx) / ady)

      ; enforce zero current since magnetic field gradient is zero at boundaries
      jxe[nx-1,*,*] = zero
      jye[*,ny-1,*] = zero
      jze[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

      ;print, 'calculating vertex-centred currents (as used by the VisIt LARE3D plugin)...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      jxv = jxe
      jyv = jye
      jzv = jze

      jxv[0:nx-2,0:ny-1,0:nz-1] = (jxe[0:nx-2,0:ny-1,0:nz-1] + jxe[1:nx-1,0:ny-1,0:nz-1]) / 2.0
      jyv[0:nx-1,0:ny-2,0:nz-1] = (jye[0:nx-1,0:ny-2,0:nz-1] + jye[0:nx-1,1:ny-1,0:nz-1]) / 2.0
      jzv[0:nx-1,0:ny-1,0:nz-2] = (jze[0:nx-1,0:ny-1,0:nz-2] + jze[0:nx-1,0:ny-1,1:nz-1]) / 2.0    
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('jx',jxv, 'jy',jyv, 'jz',jzv))
    ENDIF
    
    
    IF ((KEYWORD_SET(current) && KEYWORD_SET(bfield) && KEYWORD_SET(lorentz)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating vertex-centred Lorentz forces...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      lxv = jyv*bzv - jzv*byv
      lyv = jzv*bxv - jxv*bzv
      lzv = jxv*byv - jyv*bxv
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('lx',lxv, 'ly',lyv, 'lz',lzv))
    ENDIF
                  

    IF ((KEYWORD_SET(current) && KEYWORD_SET(resistivity)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating vertex-centred resistivities...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; initialise the 3D resistivity array
      etav = jxv
      etav[*,*,*] = eta0
      etav_anom_i = WHERE(SQRT(jxv^2 + jyv^2 + jzv^2) GE ABS(jcrit), cnt)
      IF cnt GT 0 THEN BEGIN
        x_i = etav_anom_i mod nx
        y_i = (etav_anom_i/nx) mod ny
        z_i = etav_anom_i/(nx*ny)      
        etav[x_i,y_i,z_i] = eta1      
      ENDIF
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('eta',etav))
    ENDIF
    

    IF ((KEYWORD_SET(current) && KEYWORD_SET(resistivity) && KEYWORD_SET(bfield) && KEYWORD_SET(efield)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating vertex-centred electric field components...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      exv = etav*jxv - (ds.vy*bzv - ds.vz*byv)
      eyv = etav*jyv - (ds.vz*bxv - ds.vx*bzv)
      ezv = etav*jzv - (ds.vx*byv - ds.vy*bxv)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('ex',exv, 'ey',eyv, 'ez',ezv))
    ENDIF


    IF ((KEYWORD_SET(current) && KEYWORD_SET(bfield) && KEYWORD_SET(alpha)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating vertex-centred alpha...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      alphav = jxv
      alphav[*,*,*] = 0.0
      bmag_nz_i = WHERE((bxv^2 + byv^2 + bzv^2) GT 0.0, cnt)
      IF cnt GT 0 THEN BEGIN
        x_i = bmag_nz_i mod nx
        y_i = (bmag_nz_i/nx) mod ny
        z_i = bmag_nz_i/(nx*ny)
        alphav[x_i,y_i,z_i] = (bxv[x_i,y_i,z_i]*jxv[x_i,y_i,z_i] + byv[x_i,y_i,z_i]*jyv[x_i,y_i,z_i] + bzv[x_i,y_i,z_i]*jzv[x_i,y_i,z_i])$
                             /(bxv[x_i,y_i,z_i]^2 + byv[x_i,y_i,z_i]^2 + bzv[x_i,y_i,z_i]^2)
      ENDIF
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;      

      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('alpha',alphav))
    ENDIF
    
    
    IF ((KEYWORD_SET(current) && KEYWORD_SET(resistivity) && KEYWORD_SET(bfield) && KEYWORD_SET(efield) && KEYWORD_SET(dtW)) || KEYWORD_SET(all)) THEN BEGIN 
      ;print, 'calculating face-centred electric field components...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      exf = exv
      eyf = eyv
      ezf = ezv

      exf[1:nx-1,1:ny-1,1:nz-1] = (exv[1:nx-1,1:ny-1,1:nz-1] + exv[1:nx-1,0:ny-2,1:nz-1] + exv[1:nx-1,1:ny-1,0:nz-2] + exv[1:nx-1,0:ny-2,0:nz-2]) / 4.0
      eyf[1:nx-1,1:ny-1,1:nz-1] = (eyv[1:nx-1,1:ny-1,1:nz-1] + eyv[0:nx-2,1:ny-1,1:nz-1] + eyv[1:nx-1,1:ny-1,0:nz-2] + eyv[0:nx-2,1:ny-1,0:nz-2]) / 4.0
      ezf[1:nx-1,1:ny-1,1:nz-1] = (ezv[1:nx-1,1:ny-1,1:nz-1] + ezv[0:nx-2,1:ny-1,1:nz-1] + ezv[1:nx-1,0:ny-2,1:nz-1] + ezv[0:nx-2,0:ny-2,1:nz-1]) / 4.0

      exf[1:nx-1,0,1:nz-1] = (exv[1:nx-1,0,1:nz-1] + exv[1:nx-1,0,0:nz-2]) / 2.0    
      exf[1:nx-1,1:ny-1,0] = (exv[1:nx-1,1:ny-1,0] + exv[1:nx-1,0:ny-2,0]) / 2.0    
      eyf[0,1:ny-1,1:nz-1] = (eyv[0,1:ny-1,1:nz-1] + eyv[0,1:ny-1,0:nz-2]) / 2.0    
      eyf[1:nx-1,1:ny-1,0] = (eyv[1:nx-1,1:ny-1,0] + eyv[0:nx-2,1:ny-1,0]) / 2.0    
      ezf[0,1:ny-1,1:nz-1] = (ezv[0,1:ny-1,1:nz-1] + ezv[0,0:ny-2,1:nz-1]) / 2.0
      ezf[1:nx-1,0,1:nz-1] = (ezv[1:nx-1,0,1:nz-1] + ezv[0:nx-2,0,1:nz-1]) / 2.0
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

      ;print, 'calculating vertex-centred curl of electric field...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      curlexv = curlexe
      curleyv = curleye
      curlezv = curleze

      curlexv[0:nx-2,0:ny-1,0:nz-1] = (curlexe[0:nx-2,0:ny-1,0:nz-1] + curlexe[1:nx-1,0:ny-1,0:nz-1]) / 2.0
      curleyv[0:nx-1,0:ny-2,0:nz-1] = (curleye[0:nx-1,0:ny-2,0:nz-1] + curleye[0:nx-1,1:ny-1,0:nz-1]) / 2.0
      curlezv[0:nx-1,0:ny-1,0:nz-2] = (curleze[0:nx-1,0:ny-1,0:nz-2] + curleze[0:nx-1,0:ny-1,1:nz-1]) / 2.0    
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

      ;print, 'calculating vertex-centred time derivative of the magnetic field strength...'
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      dtWvx = -(bxv*curlexv)
      dtWvy = -(byv*curleyv)
      dtWvz = -(bzv*curlezv)
      ;dtWv = -(bxv*curlexv + byv*curleyv + bzv*curlezv)

      ; enforce zero values since electric field gradient is zero at boundaries
      dtWvx[nx-1,*,*] = zero
      dtWvx[*,ny-1,*] = zero
      dtWvx[*,*,nz-1] = zero
      dtWvy[nx-1,*,*] = zero
      dtWvy[*,ny-1,*] = zero
      dtWvy[*,*,nz-1] = zero
      dtWvz[nx-1,*,*] = zero
      dtWvz[*,ny-1,*] = zero
      dtWvz[*,*,nz-1] = zero
      ;dtWv[nx-1,*,*] = zero
      ;dtWv[*,ny-1,*] = zero
      ;dtWv[*,*,nz-1] = zero
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;       

      ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('dtWx',dtWvx,'dtWy',dtWvy,'dtWz',dtWvz))
      ;ext_data = CREATE_STRUCT(ext_data, CREATE_STRUCT('dtW',dtWv))
    ENDIF        

  ENDIF ELSE BEGIN
    print, "Invalid variable passed"
    print, "Use: addcurrents, <data structure>"
    ext_data = 0
  ENDELSE  

  RETURN, ext_data

END
