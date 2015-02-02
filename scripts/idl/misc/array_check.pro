

; these variables are accessible to every function in this script
COMMON globals, nx, ny, nz



FUNCTION uxkx_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uxkx_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uxkx_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uxkx_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END


FUNCTION uxky_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uxky_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uxky_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uxky_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION uxkz_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uxkz_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uxkz_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uxkz_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION uykx_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uykx_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uykx_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uykx_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION uyky_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uyky_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uyky_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uyky_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION uykz_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uykz_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uykz_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uykz_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION uzkx_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uzkx_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uzkx_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uzkx_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION uzky_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uzky_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uzky_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uzky_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION uzkz_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+1 THEN BEGIN
    PRINT, id, ": uzkz_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+1 THEN BEGIN
    PRINT, id, ": uzkz_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+1 THEN BEGIN
    PRINT, id, ": uzkz_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END


FUNCTION dxb_chk, id, ix
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": dxb_chk failed: ix = ", ix, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION dyb_chk, id, iy
  COMMON globals

  err = 0

  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": dyb_chk failed: iy = ", iy, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION dzb_chk, id, iz
  COMMON globals

  err = 0

  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": dzb_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION dxc_chk, id, ix
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": dxc_chk failed: ix = ", ix, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION dyc_chk, id, iy
  COMMON globals

  err = 0

  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": dyc_chk failed: iy = ", iy, "."
    err = 1
  ENDIF

  RETURN, err
END

FUNCTION dzc_chk, id, iz
  COMMON globals

  err = 0

  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": dzc_chk failed: iz = ", iz, "."
    err = 1
  ENDIF

  RETURN, err
END


FUNCTION rho_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": rho_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": rho_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": rho_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION energy_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": energy_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": energy_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": energy_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION energy0_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": energy0_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": energy0_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": energy0_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION temperature_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": temperature_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": temperature_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": temperature_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION temperature0_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": temperature0_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": temperature0_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": temperature0_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END


FUNCTION radiation_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": radiation_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": radiation_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": radiation_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION alpha_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": alpha_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": alpha_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": alpha_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END

FUNCTION heat_in_chk, id, ix, iy, iz
  COMMON globals

  err = 0

  IF ix LT -1 || ix GT nx+2 THEN BEGIN
    PRINT, id, ": heat_in_chk failed: ix = ", ix, "."
    err = 1
  ENDIF
  IF iy LT -1 || iy GT ny+2 THEN BEGIN
    PRINT, id, ": heat_in_chk failed: iy = ", iy, "."
    err = 1
  ENDIF
  IF iz LT -1 || iz GT nz+2 THEN BEGIN
    PRINT, id, ": heat_in_chk failed: iz = ", iz, "."
    err = 1
  ENDIF
 
  RETURN, err
END



PRO array_check
  COMMON globals

  nx = 32
  ny = 32
  nz = 64

  FOR it=1,1,1 DO BEGIN
  
    print,"it = ",it

    z1 = 1   

    FOR redblack=1,2,1 DO BEGIN

      y1 = z1 
      
      FOR iz=1,nz,1 DO BEGIN
        
        x1 = z1
       
        FOR iy=1,ny,1 DO BEGIN

          FOR ix=x1,nx,2 DO BEGIN
              
            err = uxkx_chk(1, ix,iy,iz)
            err = dxc_chk(2, ix)
            err = dxb_chk(3, ix)
            err = uxkx_chk(4, ix-1,iy,iz)
            err = dxc_chk(5, ix-1)
            err = dxb_chk(6, ix)
            err = uyky_chk(7, ix,iy,iz)
            err = dyc_chk(8, iy)
            err = dyb_chk(9, iy)
            err = uyky_chk(10, ix,iy-1,iz)
            err = dyc_chk(11, iy-1)
            err = dyb_chk(12, iy)
            err = uzkz_chk(13, ix,iy,iz)
            err = dzc_chk(14, iz)
            err = dzb_chk(15, iz)
            err = uzkz_chk(16, ix,iy,iz-1)
            err = dzc_chk(17, iz-1)
            err = dzb_chk(18, iz)


            err = uxkx_chk(1, ix,iy,iz)
            err = temperature_chk(2, ix+1,iy,iz)
            err = dxc_chk(3, ix)
            err = dxb_chk(4, ix)
            err = uxkx_chk(5, ix-1,iy,iz)
            err = temperature_chk(6, ix-1,iy,iz)
            err = dxc_chk(7, ix-1)
            err = dxb_chk(8, ix) 
            err = uyky_chk(9, ix,iy,iz)
            err = temperature_chk(10, ix,iy+1,iz)
            err = dyc_chk(11, iy)
            err = dyb_chk(12, iy)
            err = uyky_chk(13, ix,iy-1,iz)
            err = temperature_chk(14, ix,iy-1,iz)
            err = dyc_chk(15, iy-1)
            err = dyb_chk(16, iy)  
            err = uzkz_chk(17, ix,iy,iz)
            err = temperature_chk(18, ix,iy,iz+1)
            err = dzc_chk(19, iz)
            err = dzb_chk(20, iz)
            err = uzkz_chk(21, ix,iy,iz-1)
            err = temperature_chk(22, ix,iy,iz-1)
            err = dzc_chk(23, iz-1)
            err = dzb_chk(24, iz) 

            err = uxky_chk(1, ix,iy,iz)
            err = temperature_chk(2, ix+1,iy+1,iz)
            err = temperature_chk(3, ix,iy+1,iz)
            err = temperature_chk(4, ix+1,iy-1,iz)
            err = temperature_chk(5, ix,iy-1,iz)
            err = dxb_chk(6, ix)
            err = dyc_chk(7, iy)
            err = dyc_chk(8, iy-1)  
            err = uxky_chk(9, ix-1,iy,iz)
            err = temperature_chk(10, ix,iy+1,iz)
            err = temperature_chk(11, ix-1,iy+1,iz)
            err = temperature_chk(12, ix,iy-1,iz)
            err = temperature_chk(13, ix-1,iy-1,iz)
            err = dxb_chk(14, ix)
            err = dyc_chk(15, iy)
            err = dyc_chk(16, iy-1)  
            err = uxkz_chk(17, ix,iy,iz)
            err = temperature_chk(18, ix+1,iy,iz+1)
            err = temperature_chk(19, ix,iy,iz+1)
            err = temperature_chk(20, ix+1,iy,iz-1)
            err = temperature_chk(21, ix,iy,iz-1)
            err = dxb_chk(22, ix)
            err = dzc_chk(23, iz)
            err = dzc_chk(24, iz-1)  
            err = uxkz_chk(25, ix-1,iy,iz)
            err = temperature_chk(26, ix,iy,iz+1)
            err = temperature_chk(27, ix-1,iy,iz+1)
            err = temperature_chk(28, ix,iy,iz-1)
            err = temperature_chk(29, ix-1,iy,iz-1)
            err = dxb_chk(30, ix)
            err = dzc_chk(31, iz)
            err = dzc_chk(32, iz-1)


            err = uykx_chk(1, ix,iy,iz)
            err = temperature_chk(2, ix+1,iy+1,iz)
            err = temperature_chk(3, ix+1,iy,iz)
            err = temperature_chk(4, ix-1,iy+1,iz)
            err = temperature_chk(5, ix-1,iy,iz)
            err = dyb_chk(6, iy)
            err = dxc_chk(7, ix)
            err = dxc_chk(8, ix-1) 
            err = uykx_chk(9, ix,iy-1,iz)
            err = temperature_chk(10, ix+1,iy,iz)
            err = temperature_chk(11, ix+1,iy-1,iz)
            err = temperature_chk(12, ix-1,iy,iz)
            err = temperature_chk(13, ix-1,iy-1,iz)
            err = dyb_chk(14, iy)
            err = dxc_chk(15, ix)
            err = dxc_chk(16, ix-1)
            err = uykz_chk(17, ix,iy,iz)
            err = temperature_chk(18, ix,iy+1,iz+1)
            err = temperature_chk(19, ix,iy,iz+1)
            err = temperature_chk(20, ix,iy+1,iz-1)
            err = temperature_chk(21, ix,iy,iz-1)
            err = dyb_chk(22, iy)
            err = dzc_chk(23, iz)
            err = dzc_chk(24, iz-1)
            err = uykz_chk(25, ix,iy-1,iz)
            err = temperature_chk(26, ix,iy,iz+1)
            err = temperature_chk(27, ix,iy-1,iz+1)
            err = temperature_chk(28, ix,iy,iz-1)
            err = temperature_chk(29, ix,iy-1,iz-1)
            err = dyb_chk(30,  iy)
            err = dzc_chk(31, iz)
            err = dzc_chk(32, iz-1)


            err = uzkx_chk(1, ix,iy,iz)
            err = temperature_chk(2, ix+1,iy,iz+1)
            err = temperature_chk(3, ix+1,iy,iz)
            err = temperature_chk(4, ix-1,iy,iz+1)
            err = temperature_chk(5, ix-1,iy,iz)
            err = dzb_chk(6, iz)
            err = dxc_chk(7, ix)
            err = dxc_chk(8, ix-1)
            err = uzkx_chk(9, ix,iy,iz-1)
            err = temperature_chk(10, ix+1,iy,iz)
            err = temperature_chk(11, ix+1,iy,iz-1)
            err = temperature_chk(12, ix-1,iy,iz)
            err = temperature_chk(13, ix-1,iy,iz-1)
            err = dzb_chk(14, iz)
            err = dxc_chk(15, ix)
            err = dxc_chk(16, ix-1)
            err = uzky_chk(17, ix,iy,iz)
            err = temperature_chk(18, ix,iy+1,iz+1)
            err = temperature_chk(19, ix,iy+1,iz)
            err = temperature_chk(20, ix,iy-1,iz+1)
            err = temperature_chk(21, ix,iy-1,iz)
            err = dzb_chk(22, iz)
            err = dyc_chk(23 ,iy)
            err = dyc_chk(24, iy-1)  
            err = uzky_chk(25, ix,iy,iz-1)
            err = temperature_chk(26, ix,iy+1,iz)
            err = temperature_chk(27, ix,iy+1,iz-1)
            err = temperature_chk(28, ix,iy-1,iz)
            err = temperature_chk(29, ix,iy-1,iz-1)
            err = dzb_chk(30, iz)
            err = dyc_chk(31, iy)
            err = dyc_chk(32, iy-1)

            err = radiation_chk(1, ix,iy,iz)
            err = alpha_chk(2, ix,iy,iz)
            err = temperature0_chk(3, ix,iy,iz)
            err = temperature_chk(4, ix,iy,iz)
            err = energy_chk(5, ix,iy,iz) 
            err = heat_in_chk(6, ix,iy,iz)
            err = rho_chk(7, ix, iy, iz)     
            err = energy0_chk(8, ix,iy,iz)
              
            
                                          
          ENDFOR

          x1 = 3 - x1

        ENDFOR 

        y1 = 3 - y1

      ENDFOR 

      z1 = 3 - z1
        
    ENDFOR
      
  ENDFOR

END
