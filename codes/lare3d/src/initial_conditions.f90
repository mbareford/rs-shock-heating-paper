MODULE initial_conditions

  USE shared_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: set_initial_conditions

CONTAINS
   
  SUBROUTINE set_initial_conditions
   
    INTEGER :: ix, iy ,iz

    REAL(num) :: rc, rb, rbx, rby, b_t, b_z
    REAL(num) :: dx, dy, v_p, v_r, v_t
    REAL(num) :: costh, sinth, coskz, sinkz
    
 
    ! initialise velocity arrays
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    
    ! initialise the energy, gravity and density arrays
    energy = iniEn/En0
    grav = 0.0_num
    rho = 1.0_num

    ! setup background field and apply kink perturbation      
    dx = 0.1_num * dxb(nx/2)
    dy = 0.1_num * dyb(ny/2)

    DO ix = -1,nx+2
      DO iy = -1,ny+2
        DO iz = -1,nz+2
          rc  = SQRT(xc(ix)*xc(ix)+yc(iy)*yc(iy))
          rb  = SQRT(xb(ix)*xb(ix)+yb(iy)*yb(iy))
          rbx = SQRT(xb(ix)*xb(ix)+yc(iy)*yc(iy))
          rby = SQRT(xc(ix)*xc(ix)+yb(iy)*yb(iy))

          ! background field has no z-dependence
          bx(ix,iy,iz) = calc_bx(lp_alp,lp_omg,rbx,yc(iy))
          by(ix,iy,iz) = calc_by(lp_alp,lp_omg,rby,xc(ix))          
          bz(ix,iy,iz) = calc_bz(lp_alp,lp_omg,rc)

          IF (perturbation) THEN 
            IF (rb .LE. sqrt(dx**2 + dy**2)) THEN
              costh = 1.0_num
              sinth = 0.0_num
            ELSE
              costh = xb(ix)/rb
              sinth = yb(iy)/rb
            ENDIF

            coskz = cos(ptb_k*zb(iz))
            sinkz = sin(ptb_k*zb(iz))

            b_z = calc_bz(lp_alp,lp_omg,rb)
            b_t = calc_btheta(lp_alp,lp_omg,rb)                                  
            
            v_r = 0.0_num
            v_p = 0.0_num
            v_t = 0.0_num
            IF (ABS(zb(iz)) .LE. ptb_lim) THEN
              !v_r =  ptb_sign*EXP(-4.0_num*(rb**4))*COS(pi*(zb(iz)/(2.0_num*ptb_lim)))*(costh*coskz + sinth*sinkz)

              !v_p = -ptb_sign*EXP(-4.0_num*(rb**4))*COS(pi*(zb(iz)/(2.0_num*ptb_lim)))*(sinth*coskz - costh*sinkz)&
              !               *(1.0_num - 16.0_num*(rb**4))*((b_t**2 + b_z**2)/(b_z + ptb_k*rb*b_t))

              v_r =  ptb_sign*EXP(-(rb**2)*(1.0_num + ((2.0_num*rb)**6)))*COS(pi*zb(iz)/(2.0_num*ptb_lim))*(costh*coskz + sinth*sinkz)              

              v_p = -ptb_sign*EXP(-4.0_num*(rb**4))*COS(pi*zb(iz)/(2.0_num*ptb_lim))*(sinth*coskz - costh*sinkz)&
                             *(1.0_num - 2.0_num*(rb**2) - 8.0_num*((2.0_num*rb)**6))*((b_t**2 + b_z**2)/(b_z + ptb_k*rb*b_t))

              v_t = (b_z/(b_z**2 + b_t**2))*v_p
            ENDIF 
            vx(ix,iy,iz) = ptb_amp*(v_r*costh - v_t*sinth)
            vy(ix,iy,iz) = ptb_amp*(v_r*sinth + v_t*costh)
            vz(ix,iy,iz) = ptb_amp*((-b_t/(b_z**2 + b_t**2))*v_p)

          ENDIF !perturbation
           
        ENDDO
      ENDDO
    ENDDO    
    
  END SUBROUTINE set_initial_conditions


  !----------------------------------------------------------------------------------
  !These functions calculate magnetic field components at specific radial coordinates
  !----------------------------------------------------------------------------------

  FUNCTION calc_bx(alp,omg,r,y)
    REAL(num), INTENT(IN) :: alp,omg,r,y
    REAL(num) :: calc_bx

    IF (r.EQ.0.0_num) THEN
      calc_bx = 0.0_num
    ELSE
      calc_bx = calc_btheta(alp,omg,r)*(-y/r)
    ENDIF
  END FUNCTION calc_bx

  FUNCTION calc_by(alp,omg,r,x)    
    REAL(num), INTENT(IN) :: alp,omg,r,x
    REAL(num) :: calc_by

    IF (r.EQ.0.0_num) THEN
      calc_by = 0.0_num
    ELSE
      calc_by = calc_btheta(alp,omg,r)*(x/r)
    ENDIF
  END FUNCTION calc_by

  FUNCTION calc_btheta(alp,omg,r)    
    REAL(num), INTENT(IN) :: alp,omg,r
    REAL(num) :: calc_btheta

    IF (r.EQ.0.0_num) THEN
      calc_btheta = 0.0_num
    ELSEIF (r.LT.lp_rad) THEN
      calc_btheta = alp*r*((1.0_num - r**2)**omg)
    ELSE
      calc_btheta = 0.0_num
    ENDIF
  END FUNCTION calc_btheta

  FUNCTION calc_bz(alp,omg,r)
    REAL(num), INTENT(IN) :: alp,omg,r
    REAL(num) :: calc_bz
   
    IF (r.EQ.0.0_num) THEN
      calc_bz = 1.0_num
    ELSEIF (r.LT.lp_rad) THEN
      calc_bz = SQRT(1.0_num&
                     + ((alp**2)/(2.0_num*omg + 1.0_num))*((1.0_num - r**2)**(2.0_num*omg + 1.0_num) - 1.0_num)&
                     - (alp**2)*(r**2)*((1.0_num - r**2)**(2.0_num*omg)))
    ELSE
      calc_bz = SQRT(1.0_num - (alp**2)/(2.0_num*omg + 1.0_num))
    END IF
  END FUNCTION calc_bz  

END MODULE initial_conditions

  
