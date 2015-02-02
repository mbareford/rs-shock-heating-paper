PRO normaliser, units

  mu0 = 4.0d-7*!dpi
  mp = 1.6726d-27
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

  P0 = B0^2/mu0    ; N m^-2
  e0 = P0/rho0     ; J 
  T0 = e0*(m0/kB)  ; K  

  ; conduction kappa
  K0 = 1.0d-11*(T0^(3.5d0))*tA/(rho0*(L0^2)*e0)

  IF STRCMP(units,'cgs') EQ 1 THEN BEGIN
    print,'Length, L0 = ', L0*1.0d2, ' cm'
    print,'Number density, n0 = ', n0/1.0d6, ' cm^-3'
    print,'Magnetic field strength, B0 = ', B0*1.0d4, ' G'
    print,''
    print,'Density, rho0 = ', rho0/(1.0d3), ' g cm^-3'
    print,'Alfven velocity, vA = ', vA*1.0d2, ' cm s^-1'
    print,'Alfven time, tA = ', tA, ' s'
    print,''
    print,'Pressure, P0 = ', P0*1.0d1, ' dyn cm^-2'
    print,'Energy, e0 = ', e0*1.0d7, ' erg'
    print,'Temperature, T0 = ', T0, ' K'
  ENDIF ELSE BEGIN
    print,'Length, L0 = ', L0, ' m'
    print,'Number density, n0 = ', n0, ' m^-3'
    print,'Magnetic field strength, B0 = ', B0, ' T'
    print,''
    print,'Density, rho0 = ', rho0, ' kg m^-3'
    print,'Alfven velocity, vA = ', vA, ' m s^-1'
    print,'Alfven time, tA = ', tA, ' s'
    print,''
    print,'Pressure, P0 = ', P0, ' N m^-2'
    print,'Energy, e0 = ', e0, ' J'
    print,'Temperature, T0 = ', T0, ' K'
  ENDELSE
 
END
