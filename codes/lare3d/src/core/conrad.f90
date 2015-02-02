MODULE conrad

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PUBLIC :: conduct, con_x_loss, con_y_loss, con_z_loss, radiate, rad_dt
  
CONTAINS

  ! calculate the energy lost due to conduction through the x boundaries
  ! this function is called by output_routines in diagnostics.F90
  FUNCTION con_x_loss
    
    REAL(num) :: con_x_loss, flux 
    REAL(num) :: bxc, byc, bzc, bpx, txb
    REAL(num) :: ux, uy, uz
    REAL(num) :: temp_fac, flux_sign, flux_dt

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: uxkx, uxky, uxkz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp

         
    ix = -1
    IF (xbc_min .EQ. BC_OTHER) THEN
      IF (proc_x_min .EQ. MPI_PROC_NULL) THEN
        ix = 0
      ELSE IF (proc_x_max .EQ. MPI_PROC_NULL) THEN
        ix = nx
      END IF
    END IF

    IF (ix .EQ. -1) THEN
      con_x_loss = 0.0_num
      RETURN
    END IF

    
    ALLOCATE(temp(-1:nx+2,-1:ny+2,-1:nz+2))
    temp_fac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) temp_fac = 2.0_num
    temp = temp_fac*(gamma - 1.0_num)*(energy - (1.0_num - xi_n)*ionise_pot)/(2.0_num - xi_n)
  
    ALLOCATE(uxkx(1:ny,1:nz), uxky(1:ny,1:nz), uxkz(1:ny,1:nz))

    uxkx = 0.0_num
    uxky = 0.0_num
    uxkz = 0.0_num  


    DO iz=1,nz
      DO iy=1,ny

        bxc = bx(ix,iy,iz) 
        byc = (by(ix,iy,iz) + by(ix+1,iy,iz) + by(ix,iy-1,iz) + by(ix+1,iy-1,iz))/4.0_num
        bzc = (bz(ix,iy,iz) + bz(ix+1,iy,iz) + bz(ix,iy,iz-1) + bz(ix+1,iy,iz-1))/4.0_num

        bpx = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
        bpx = MAX(bpx, none_zero) 

        ux = bxc/bpx       
        uy = byc/bpx       
        uz = bzc/bpx       

        txb = (temp(ix,iy,iz) + temp(ix+1,iy,iz))/2.0_num

        uxkx(iy,iz) = ux*ux*K0*(txb**2.5_num)
        uxky(iy,iz) = ux*uy*K0*(txb**2.5_num)
        uxkz(iy,iz) = ux*uz*K0*(txb**2.5_num) 

        uxkx(iy,iz) = uxkx(iy,iz) + ((con_b_min**2)*K0*(txb**2.5_num))/(bpx**2 + con_b_min**2) 

      END DO
    END DO
    
       
    flux_sign = 1.0_num
    IF (ix .EQ. nx) flux_sign = -1.0_num

    flux_dt = dt
        

    con_x_loss = 0.0_num

    DO iz=1,nz
      DO iy=1,ny

        flux = + uxkx(iy,iz)*(temp(ix+1,iy,iz) - temp(ix,iy,iz))/dxc(ix) &
               + uxky(iy,iz)*0.5_num*(temp(ix+1,iy+1,iz) - temp(ix+1,iy-1,iz) + temp(ix,iy+1,iz) - temp(ix,iy-1,iz))/(dyc(iy-1)+dyc(iy)) &
               + uxkz(iy,iz)*0.5_num*(temp(ix+1,iy,iz+1) - temp(ix+1,iy,iz-1) + temp(ix,iy,iz+1) - temp(ix,iy,iz-1))/(dzc(iz-1)+dzc(iz))
    
        con_x_loss = con_x_loss + flux_sign*flux*flux_dt*(dyb(iy)*dzb(iz))
        
      END DO
    END DO    


    DEALLOCATE(temp)
    DEALLOCATE(uxkx,uxky,uxkz)
  
  END FUNCTION con_x_loss

  ! calculate the energy lost due to conduction through the y boundaries
  ! this function is called by output_routines in diagnostics.F90
  FUNCTION con_y_loss
    
    REAL(num) :: con_y_loss, flux 
    REAL(num) :: bxc, byc, bzc, bpy, tyb
    REAL(num) :: ux, uy, uz
    REAL(num) :: temp_fac, flux_sign, flux_dt

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: uykx, uyky, uykz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp
    
         
    iy = -1
    IF (ybc_min .EQ. BC_OTHER) THEN
      IF (proc_y_min .EQ. MPI_PROC_NULL) THEN
        iy = 0
      ELSE IF (proc_y_max .EQ. MPI_PROC_NULL) THEN
        iy = ny
      END IF
    END IF

    IF (iy .EQ. -1) THEN
      con_y_loss = 0.0_num
      RETURN
    END IF

    
    ALLOCATE(temp(-1:nx+2,-1:ny+2,-1:nz+2))
    temp_fac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) temp_fac = 2.0_num
    temp = temp_fac*(gamma - 1.0_num)*(energy - (1.0_num - xi_n)*ionise_pot)/(2.0_num - xi_n)
  
    ALLOCATE(uykx(1:nx,1:nz), uyky(1:nx,1:nz), uykz(1:nx,1:nz))
    
    uykx = 0.0_num
    uyky = 0.0_num
    uykz = 0.0_num  
    

    DO iz=1,nz
      DO ix=1,nx

        bxc = (bx(ix,iy,iz) + bx(ix,iy+1,iz) + bx(ix-1,iy,iz) + bx(ix-1,iy+1,iz))/4.0_num
        byc = by(ix,iy,iz)      
        bzc = (bz(ix,iy,iz) + bz(ix,iy+1,iz) + bz(ix,iy,iz-1) + bz(ix,iy+1,iz-1))/4.0_num
          
        bpy = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
        bpy = MAX(bpy, none_zero) 
          
        ux = bxc/bpy      
        uy = byc/bpy       
        uz = bzc/bpy      
          
        tyb = (temp(ix,iy,iz) + temp(ix,iy+1,iz))/2.0_num

        uykx(ix,iz) = uy*ux*K0*(tyb**2.5_num)
        uyky(ix,iz) = uy*uy*K0*(tyb**2.5_num)  
        uykz(ix,iz) = uy*uz*K0*(tyb**2.5_num)   

        uyky(ix,iz) = uyky(ix,iz) + ((con_b_min**2)*K0*(tyb**2.5_num))/(bpy**2 + con_b_min**2)  

      END DO
    END DO

        
    flux_sign = 1.0_num
    IF (iy .EQ. ny) flux_sign = -1.0_num

    flux_dt = dt
    

    con_y_loss = 0.0_num    

    DO iz=1,nz
      DO ix=1,nx

        flux =   uykx(ix,iz)*0.5_num*(temp(ix+1,iy+1,iz) - temp(ix-1,iy+1,iz) + temp(ix+1,iy,iz) - temp(ix-1,iy,iz))/(dxc(ix-1)+dxc(ix)) &
               + uyky(ix,iz)*(temp(ix,iy+1,iz) - temp(ix,iy,iz))/dyc(iy) &
               + uykz(ix,iz)*0.5_num*(temp(ix,iy+1,iz+1) - temp(ix,iy+1,iz-1) + temp(ix,iy,iz+1) - temp(ix,iy,iz-1))/(dzc(iz-1)+dzc(iz))

        con_y_loss = con_y_loss + flux_sign*flux*flux_dt*(dxb(ix)*dzb(iz))
        
      END DO
    END DO    


    DEALLOCATE(temp)
    DEALLOCATE(uykx,uyky,uykz)
    
  END FUNCTION con_y_loss

  ! calculate the energy lost due to conduction through the z boundaries
  ! this function is called by output_routines in diagnostics.F90
  FUNCTION con_z_loss
    
    REAL(num) :: con_z_loss, flux 
    REAL(num) :: bxc, byc, bzc, bpz, tzb
    REAL(num) :: ux, uy, uz
    REAL(num) :: temp_fac, flux_sign, flux_dt

    REAL(num), DIMENSION(:, :), ALLOCATABLE :: uzkx, uzky, uzkz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: temp

             
    iz = -1
    IF (zbc_min .EQ. BC_OTHER) THEN
      IF (proc_z_min .EQ. MPI_PROC_NULL) THEN
        iz = 0
      ELSE IF (proc_z_max .EQ. MPI_PROC_NULL) THEN
        iz = nz
      END IF
    END IF

    IF (iz .EQ. -1) THEN
      con_z_loss = 0.0_num
      RETURN
    END IF

    
    ALLOCATE(temp(-1:nx+2,-1:ny+2,-1:nz+2))
    temp_fac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) temp_fac = 2.0_num
    temp = temp_fac*(gamma - 1.0_num)*(energy - (1.0_num - xi_n)*ionise_pot)/(2.0_num - xi_n)
  
    ALLOCATE(uzkx(1:nx,1:ny), uzky(1:nx,1:ny), uzkz(1:nx,1:ny))

    uzkx = 0.0_num
    uzky = 0.0_num
    uzkz = 0.0_num  


    DO iy=1,ny
      DO ix=1,nx

        bxc = (bx(ix,iy,iz) + bx(ix,iy,iz+1) + bx(ix-1,iy,iz) + bx(ix-1,iy,iz+1))/4.0_num
        byc = (by(ix,iy,iz) + by(ix,iy,iz+1) + by(ix,iy-1,iz) + by(ix,iy-1,iz+1))/4.0_num
        bzc = bz(ix,iy,iz)      

        bpz = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
        bpz = MAX(bpz, none_zero) 
          
        ux = bxc/bpz
        uy = byc/bpz
        uz = bzc/bpz

        tzb = (temp(ix,iy,iz) + temp(ix,iy,iz+1))/2.0_num
        
        uzkx(ix,iy) = uz*ux*K0*(tzb**(2.5_num))
        uzky(ix,iy) = uz*uy*K0*(tzb**(2.5_num))
        uzkz(ix,iy) = uz*uz*K0*(tzb**(2.5_num))

        uzkz(ix,iy) = uzkz(ix,iy) + ((con_b_min**2)*K0*(tzb**2.5_num))/(bpz**2 + con_b_min**2)

      END DO
    END DO

    
    flux_sign = 1.0_num
    IF (iz .EQ. nz) flux_sign = -1.0_num

    flux_dt = dt
            

    con_z_loss = 0.0_num

    DO iy=1,ny
      DO ix=1,nx

        flux =   uzkx(ix,iy)*0.5_num*(temp(ix+1,iy,iz+1) - temp(ix-1,iy,iz+1) + temp(ix+1,iy,iz) - temp(ix-1,iy,iz))/(dxc(ix-1)+dxc(ix)) &
               + uzky(ix,iy)*0.5_num*(temp(ix,iy+1,iz+1) - temp(ix,iy-1,iz+1) + temp(ix,iy+1,iz) - temp(ix,iy-1,iz))/(dyc(iy-1)+dyc(iy)) &        
               + uzkz(ix,iy)*(temp(ix,iy,iz+1) - temp(ix,iy,iz))/dzc(iz)

        con_z_loss = con_z_loss + flux_sign*flux*flux_dt*(dxb(ix)*dyb(iy))
        
      END DO
    END DO    


    DEALLOCATE(temp)
    DEALLOCATE(uzkx,uzky,uzkz)    
  
  END FUNCTION con_z_loss




  SUBROUTINE conduct

    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uxkx, uxky, uxkz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uykx, uyky, uykz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: uzkx, uzky, uzkz
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: limit, rad_a1, rad_a2
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: energy0, temp
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: local_energy_supp

    REAL(num) :: txb, tyb, tzb
    REAL(num) :: bxc, byc, bzc, bpx, bpy, bpz
    REAL(num) :: ux, uy, uz
    REAL(num) :: q_shx, q_shy, q_shz, q_sh, q_f, q_nl
    REAL(num) :: Lr, Lr_enbar, alpha, rad_en
    REAL(num) :: a1, a2, max_energy0, glb_max_energy0, temp_fac
    REAL(num) :: rsd, it_err, max_err, glb_max_err, glb_max_it_cnt
    REAL(num) :: it_cnt_rl, loss_rad0, energymin
    REAL(num) :: energyn, rad_min_en

    INTEGER :: redblack, x1, z1, it_cnt

    LOGICAL :: converged
    
    
    ALLOCATE(uxkx(-1:nx+1,-1:ny+1,-1:nz+1), uxky(-1:nx+1,-1:ny+1,-1:nz+1), uxkz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uykx(-1:nx+1,-1:ny+1,-1:nz+1), uyky(-1:nx+1,-1:ny+1,-1:nz+1), uykz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(uzkx(-1:nx+1,-1:ny+1,-1:nz+1), uzky(-1:nx+1,-1:ny+1,-1:nz+1), uzkz(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(rad_a1(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(rad_a2(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(energy0(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(temp(-1:nx+2,-1:ny+2,-1:nz+2))
    ALLOCATE(local_energy_supp(1:nx,1:ny,1:nz))
    
    uxkx = 0.0_num
    uxky = 0.0_num
    uxkz = 0.0_num
    uykx = 0.0_num
    uyky = 0.0_num
    uykz = 0.0_num
    uzkx = 0.0_num
    uzky = 0.0_num
    uzkz = 0.0_num

    max_energy0 = MAXVAL(energy)
    CALL MPI_ALLREDUCE(max_energy0, glb_max_energy0, 1, mpireal, MPI_MAX, comm, errcode)
    energy0 = energy

    temp_fac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) temp_fac = 2.0_num
    temp = temp_fac*(gamma - 1.0_num)*(energy - (1.0_num - xi_n)*ionise_pot)/(2.0_num - xi_n)
    energymin = MINVAL(energy)
                     
    IF (conduction) THEN
      
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1 
  
            ! x-face-centred b field
            bxc = bx(ix,iy,iz) 
            byc = (by(ix,iy,iz) + by(ix+1,iy,iz) + by(ix,iy-1,iz) + by(ix+1,iy-1,iz))/4.0_num
            bzc = (bz(ix,iy,iz) + bz(ix+1,iy,iz) + bz(ix,iy,iz-1) + bz(ix+1,iy,iz-1))/4.0_num
            bpx = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
            bpx = MAX(bpx, none_zero) 
            ! temperature at x boundary
            txb = (temp(ix,iy,iz) + temp(ix+1,iy,iz))/2.0_num
            ! direction of magnetic field on x face 
            ux = bxc/bpx       
            uy = byc/bpx       
            uz = bzc/bpx       
            ! kappa along magnetic field, now a vector  
            uxkx(ix,iy,iz) = ux*ux*K0*(txb**2.5_num)
            uxky(ix,iy,iz) = ux*uy*K0*(txb**2.5_num)
            uxkz(ix,iy,iz) = ux*uz*K0*(txb**2.5_num) 
            ! add symmetric conduction near b=0 points 
            uxkx(ix,iy,iz) = uxkx(ix,iy,iz) + ((con_b_min**2)*K0*(txb**2.5_num))/(bpx**2 + con_b_min**2) 
  
            ! y-face-centred b field
            bxc = (bx(ix,iy,iz) + bx(ix,iy+1,iz) + bx(ix-1,iy,iz) + bx(ix-1,iy+1,iz))/4.0_num
            byc = by(ix,iy,iz)      
            bzc = (bz(ix,iy,iz) + bz(ix,iy+1,iz) + bz(ix,iy,iz-1) + bz(ix,iy+1,iz-1))/4.0_num
            bpy = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
            bpy = MAX(bpy, none_zero) 
            ! temperature at y boundary
            tyb = (temp(ix,iy,iz) + temp(ix,iy+1,iz))/2.0_num
            ! direction of magnetic field on y face 
            ux = bxc/bpy      
            uy = byc/bpy       
            uz = bzc/bpy      
            ! kappa along magnetic field, now a vector  
            uykx(ix,iy,iz) = uy*ux*K0*(tyb**2.5_num)
            uyky(ix,iy,iz) = uy*uy*K0*(tyb**2.5_num)  
            uykz(ix,iy,iz) = uy*uz*K0*(tyb**2.5_num)   
            ! add symmetic conduction near b=0 points 
            uyky(ix,iy,iz) = uyky(ix,iy,iz) + ((con_b_min**2)*K0*(tyb**2.5_num))/(bpy**2 + con_b_min**2)  
  
            ! z-face-centred B field
            bxc = (bx(ix,iy,iz) + bx(ix,iy,iz+1) + bx(ix-1,iy,iz) + bx(ix-1,iy,iz+1))/4.0_num
            byc = (by(ix,iy,iz) + by(ix,iy,iz+1) + by(ix,iy-1,iz) + by(ix,iy-1,iz+1))/4.0_num         
            bzc = bz(ix,iy,iz) 
            bpz = SQRT(bxc**2 + byc**2 + bzc**2 + con_b_min**2)
            bpz = MAX(bpz, none_zero) 
            ! temperature at z boundary
            tzb = (temp(ix,iy,iz) + temp(ix,iy,iz+1))/2.0_num 
            ! direction of magnetic field on z face 
            ux = bxc/bpz      
            uy = byc/bpz       
            uz = bzc/bpz      
            ! kappa along magnetic field, now a vector  
            uzkx(ix,iy,iz) = uz*ux*K0*(tzb**2.5_num)
            uzky(ix,iy,iz) = uz*uy*K0*(tzb**2.5_num)   
            uzkz(ix,iy,iz) = uz*uz*K0*(tzb**2.5_num)   
            ! add symmetic conduction near b=0 points 
            uzkz(ix,iy,iz) = uzkz(ix,iy,iz) + ((con_b_min**2)*K0*(tzb**2.5_num))/(bpz**2 + con_b_min**2)   
  
          END DO
        END DO  
      END DO
    END IF ! end of if conduction clause


    IF (conduction .AND. con_flux_limited) THEN
      ALLOCATE(limit(-1:nx+2,-1:ny+2,-1:nz+2))

      DO iz = 0, nz + 1
        DO iy = 0, ny + 1
          DO ix = 0, nx + 1  
            ! estimate the parallel heat flux at the centre of a cell
            q_shx =   (uxkx(ix,iy,iz) + uxkx(ix-1,iy,iz))*(temp(ix+1,iy,iz) - temp(ix-1,iy,iz))/dxc(ix) & 
                    + (uxky(ix,iy,iz) + uxky(ix,iy-1,iz))*(temp(ix,iy+1,iz) - temp(ix,iy-1,iz))/dyc(iy) &
                    + (uxkz(ix,iy,iz) + uxkz(ix,iy,iz-1))*(temp(ix,iy,iz+1) - temp(ix,iy,iz-1))/dzc(iz) 

            q_shy =   (uykx(ix,iy,iz) + uykx(ix-1,iy,iz))*(temp(ix+1,iy,iz) - temp(ix-1,iy,iz))/dxc(ix) & 
                    + (uyky(ix,iy,iz) + uyky(ix,iy-1,iz))*(temp(ix,iy+1,iz) - temp(ix,iy-1,iz))/dyc(iy) &
                    + (uykz(ix,iy,iz) + uykz(ix,iy,iz-1))*(temp(ix,iy,iz+1) - temp(ix,iy,iz-1))/dzc(iz) 

            q_shz =   (uzkx(ix,iy,iz) + uzkx(ix-1,iy,iz))*(temp(ix+1,iy,iz) - temp(ix-1,iy,iz))/dxc(ix) & 
                    + (uzky(ix,iy,iz) + uzky(ix,iy-1,iz))*(temp(ix,iy+1,iz) - temp(ix,iy-1,iz))/dyc(iy) &
                    + (uzkz(ix,iy,iz) + uzkz(ix,iy,iz-1))*(temp(ix,iy,iz+1) - temp(ix,iy,iz-1))/dzc(iz) 

            q_sh = SQRT((q_shx**2 + q_shy**2 + q_shz**2)/16.0_num)

            ! estimate the free streaming limit
            q_f = SQRT(mp_si/me_si)*con_limit*rho(ix,iy,iz)*(MIN(temp(ix,iy,iz), 1.e8_num/Tp0)**1.5_num)
            IF ((eos_number .EQ. EOS_IDEAL) .AND. (.NOT. neutral_gas)) THEN
              ! use the electron number density
              q_f = q_f / 2.0_num  
            END IF

            q_nl = 1.0_num / (1.0_num/MAX(q_sh, none_zero) + 1.0_num/MAX(q_f, none_zero))

            limit(ix,iy,iz) = q_nl/(2.0_num*MAX(q_sh, none_zero))
          END DO
        END DO  
      END DO  

      DO iz = 0, nz+1
        DO iy = 0, ny+1
          DO ix = 0, nx+1  
            uxkx(ix,iy,iz) = uxkx(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix+1,iy,iz)) 
            uxky(ix,iy,iz) = uxky(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix+1,iy,iz))
            uxkz(ix,iy,iz) = uxkz(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix+1,iy,iz))
            uykx(ix,iy,iz) = uykx(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy+1,iz))
            uyky(ix,iy,iz) = uyky(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy+1,iz)) 
            uykz(ix,iy,iz) = uykz(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy+1,iz))
            uzkx(ix,iy,iz) = uzkx(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy,iz+1)) 
            uzky(ix,iy,iz) = uzky(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy,iz+1))
            uzkz(ix,iy,iz) = uzkz(ix,iy,iz)*(limit(ix,iy,iz) + limit(ix,iy,iz+1))
          END DO
        END DO
      END DO

      DEALLOCATE(limit)
    END IF ! end of if conduction .AND. con_flux_limited clause
       
    
    rad_a1 = 0.0_num
    rad_a2 = 0.0_num
    IF (radiation .AND. rad_implicit) THEN
      loss_rad0 = loss_rad
      rad_min_en = ((2.0_num/(gamma-1.0_num))*((kb_si*rad_min_tp)/m0))/En0

      ! calculate the radiation losses and the alpha exponent      
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx

            CALL calc_rad_loss(rho(ix,iy,iz), energy(ix,iy,iz), xi_n(ix,iy,iz), Lr, alpha)
            
            ! express Lr in terms of normalised energy units
            Lr_enbar = Lr*((dt*tm0)/(rho(ix,iy,iz)*Dn0))/En0
            rad_a1(ix,iy,iz) = alpha*Lr_enbar/energy(ix,iy,iz)
            rad_a2(ix,iy,iz) = Lr_enbar*(alpha - 1.0_num)

            rad_en = energy(ix,iy,iz) - (energy(ix,iy,iz) + rad_a2(ix,iy,iz))/(1.0_num + rad_a1(ix,iy,iz))
            IF (rad_en .GT. 0.0_num) THEN
              IF (rad_min_en .GT. (energy(ix,iy,ix)-rad_en)) THEN
                rad_en = energy(ix,iy,iz) - rad_min_en
              END IF

              loss_rad = loss_rad + rad_en*rho(ix,iy,iz)*cv(ix,iy,iz)
            END IF
  
          END DO
        END DO    		
      END DO
    END IF


    converged = .FALSE.
    local_energy_supp = 0.0_num
    
    ! iterate to get energy^{n+1} by SOR Gauss-Seidel 		
    iterate: DO it_cnt = 1, sor_it_lim
      max_err = 0.0_num
      it_err = 0.0_num      

      z1 = 1 
      DO redblack = 1, 2

        DO iz = 1, nz
          x1 = z1
          DO iy = 1, ny

             DO ix = x1, nx, 2
              
               a1 = 0.0_num
               a2 = 0.0_num

               ! d^2/dx^2, d^2/dy^2 and d^2/dz^2 derivatives, without temperature
               a1 = a1 + uxkx(ix,iy,iz)/(dxc(ix)*dxb(ix)) + uxkx(ix-1,iy,iz)/(dxc(ix-1)*dxb(ix))
               a1 = a1 + uyky(ix,iy,iz)/(dyc(iy)*dyb(iy)) + uyky(ix,iy-1,iz)/(dyc(iy-1)*dyb(iy))
               a1 = a1 + uzkz(ix,iy,iz)/(dzc(iz)*dzb(iz)) + uzkz(ix,iy,iz-1)/(dzc(iz-1)*dzb(iz))  
              
               ! d^2/dx^2, d^2/dy^2 and d^2/dz^2 derivatives
               a2 = a2 + uxkx(ix,iy,iz)*temp(ix+1,iy,iz)/(dxc(ix)*dxb(ix)) + uxkx(ix-1,iy,iz)*temp(ix-1,iy,iz)/(dxc(ix-1)*dxb(ix)) 
               a2 = a2 + uyky(ix,iy,iz)*temp(ix,iy+1,iz)/(dyc(iy)*dyb(iy)) + uyky(ix,iy-1,iz)*temp(ix,iy-1,iz)/(dyc(iy-1)*dyb(iy))              
               a2 = a2 + uzkz(ix,iy,iz)*temp(ix,iy,iz+1)/(dzc(iz)*dzb(iz)) + uzkz(ix,iy,iz-1)*temp(ix,iy,iz-1)/(dzc(iz-1)*dzb(iz))               
                
               ! d^2/dxdy cross derivatives                  
               a2 = a2 + uxky(ix,iy,iz)*(temp(ix+1,iy+1,iz) + temp(ix,iy+1,iz) - temp(ix+1,iy-1,iz) - temp(ix,iy-1,iz)) / (2.0_num*dxb(ix)*(dyc(iy) + dyc(iy-1)))  
               a2 = a2 - uxky(ix-1,iy,iz)*(temp(ix,iy+1,iz) + temp(ix-1,iy+1,iz) - temp(ix,iy-1,iz) - temp(ix-1,iy-1,iz)) / (2.0_num*dxb(ix)*(dyc(iy) + dyc(iy-1)))  

               ! d^2/dxdz cross derivatives                  
               a2 = a2 + uxkz(ix,iy,iz)*(temp(ix+1,iy,iz+1) + temp(ix,iy,iz+1) - temp(ix+1,iy,iz-1) - temp(ix,iy,iz-1)) / (2.0_num*dxb(ix)*(dzc(iz) + dzc(iz-1)))  
               a2 = a2 - uxkz(ix-1,iy,iz)*(temp(ix,iy,iz+1) + temp(ix-1,iy,iz+1) - temp(ix,iy,iz-1) - temp(ix-1,iy,iz-1)) / (2.0_num*dxb(ix)*(dzc(iz) + dzc(iz-1)))
        
               ! d^2/dydx cross derivatives
               a2 = a2 + uykx(ix,iy,iz)*(temp(ix+1,iy+1,iz) + temp(ix+1,iy,iz) - temp(ix-1,iy+1,iz) - temp(ix-1,iy,iz)) / (2.0_num*dyb(iy)*(dxc(ix) + dxc(ix-1)))  
               a2 = a2 - uykx(ix,iy-1,iz)*(temp(ix+1,iy,iz) + temp(ix+1,iy-1,iz) - temp(ix-1,iy,iz) - temp(ix-1,iy-1,iz)) / (2.0_num*dyb(iy)*(dxc(ix) + dxc(ix-1)))   
                  
               ! d^2/dydz cross derivatives
               a2 = a2 + uykz(ix,iy,iz)*(temp(ix,iy+1,iz+1) + temp(ix,iy,iz+1) - temp(ix,iy+1,iz-1) - temp(ix,iy,iz-1)) / (2.0_num*dyb(iy)*(dzc(iz) + dzc(iz-1)))  
               a2 = a2 - uykz(ix,iy-1,iz)*(temp(ix,iy,iz+1) + temp(ix,iy-1,iz+1) - temp(ix,iy,iz-1) - temp(ix,iy-1,iz-1)) / (2.0_num*dyb(iy)*(dzc(iz) + dzc(iz-1)))
              
               ! d^2/dzdx cross derivatives
               a2 = a2 + uzkx(ix,iy,iz)*(temp(ix+1,iy,iz+1) + temp(ix+1,iy,iz) - temp(ix-1,iy,iz+1) - temp(ix-1,iy,iz)) / (2.0_num*dzb(iz)*(dxc(ix) + dxc(ix-1)))  
               a2 = a2 - uzkx(ix,iy,iz-1)*(temp(ix+1,iy,iz) + temp(ix+1,iy,iz-1) - temp(ix-1,iy,iz) - temp(ix-1,iy,iz-1)) / (2.0_num*dzb(iz)*(dxc(ix) + dxc(ix-1)))  

               ! d^2/dzdy cross derivatives
               a2 = a2 + uzky(ix,iy,iz)*(temp(ix,iy+1,iz+1) + temp(ix,iy+1,iz) - temp(ix,iy-1,iz+1) - temp(ix,iy-1,iz)) / (2.0_num*dzb(iz)*(dyc(iy) + dyc(iy-1)))  
               a2 = a2 - uzky(ix,iy,iz-1)*(temp(ix,iy+1,iz) + temp(ix,iy+1,iz-1) - temp(ix,iy-1,iz) - temp(ix,iy-1,iz-1)) / (2.0_num*dzb(iz)*(dyc(iy) + dyc(iy-1)))     

                              
               a1 = (dt/rho(ix, iy, iz))*(temp(ix,iy,iz)/energy(ix,iy,iz))*a1 + rad_a1(ix,iy,iz)
               a2 = (dt/rho(ix, iy, iz))*a2 + rad_a2(ix,iy,iz)  
                 
               rsd = energy(ix,iy,iz) - (energy0(ix,iy,iz) + a2)/(1.0_num + a1)

               !energy(ix,iy,iz) = MAX(energy(ix,iy,iz) - sor_w*rsd, 0.0_num)
               energyn = energy(ix,iy,iz) - sor_w*rsd
               IF (energyn .LT. energymin) THEN
                 ! record the energy supplement required to prevent energies falling below minimum value                 
                 local_energy_supp(ix,iy,iz) = local_energy_supp(ix,iy,iz) + (energymin - energyn)
                 energy(ix,iy,iz) = energymin
               ELSE
                 energy(ix,iy,iz) = energyn
               END IF
                 
               temp(ix,iy,iz) = temp_fac*(gamma - 1.0_num)*(energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz))*ionise_pot)/(2.0_num - xi_n(ix,iy,iz))
                
               IF (sor_err_wgt == 1) THEN
                 it_err = ABS(rsd)/glb_max_energy0
               ELSE IF (sor_err_wgt == 2) THEN
                 it_err = ABS(rsd)/energy0(ix,iy,iz) 
               ELSE IF (sor_err_wgt == 3) THEN
                 it_err = ABS(rsd)/energy(ix,iy,iz)
               ELSE
                 it_err = ABS(rsd)
               END IF
               max_err = MAX(max_err, it_err)               
               
            END DO ! end of (DO ix = x1, nx, 2) loop
            x1 = 3 - x1
          END DO ! end of (DO iy = 1, ny) loop 
          
        END DO ! end of (DO iz = 1, nz) loop 
        z1 = 3 - z1
        
        CALL energy_bcs
        temp = temp_fac*(gamma - 1.0_num)*(energy - (1.0_num - xi_n)*ionise_pot)/(2.0_num - xi_n)

      END DO ! end of (DO redblack = 1, 2) loop
      
      it_cnt_rl = REAL(it_cnt)
      CALL MPI_ALLREDUCE(it_cnt_rl, glb_max_it_cnt, 1, mpireal, MPI_MAX, comm, errcode)
      CALL MPI_ALLREDUCE(max_err, glb_max_err, 1, mpireal, MPI_MAX, comm, errcode)
            
      IF (glb_max_err .LT. sor_err) THEN
        converged = .TRUE.  
        EXIT iterate
      END IF
           
    END DO iterate
    
    IF (rank == 0) THEN
      PRINT *, time, " ", INT(glb_max_it_cnt), " ", glb_max_err
    END IF

    IF (converged) THEN
      con_supp = con_supp + SUM(local_energy_supp(1:nx,1:ny,1:nz)*rho(1:nx,1:ny,1:nz)*cv(1:nx,1:ny,1:nz))       
    ELSE
      IF (rank == 0) THEN          
        PRINT * , "Conduction and/or implicit radiation failed at t = ", time, "."
      END IF
      ! restore the energy and turn off conduction and implicit radiation
      energy = energy0
      conduction = .FALSE.
      rad_implicit = .FALSE.      
    END IF

    DEALLOCATE(uxkx, uxky, uxkz)
    DEALLOCATE(uykx, uyky, uykz)
    DEALLOCATE(uzkx, uzky, uzkz)
    DEALLOCATE(rad_a1, rad_a2)
    DEALLOCATE(energy0)
    DEALLOCATE(temp)
    DEALLOCATE(local_energy_supp)
    
  END SUBROUTINE conduct


  ! calculate the radiation loss according to the method given by                               
  ! Rosner, R., Tucker, W. H., & Vaiana, G. S., 1978a, ApJ, 220, 643
  SUBROUTINE calc_rtv_rad_loss(rhobar, enbar, nf, rad_loss, rad_alpha)  
    REAL(num), INTENT(IN) :: rhobar, enbar, nf
    REAL(num), INTENT(OUT) :: rad_loss, rad_alpha 

    REAL(num), PARAMETER :: n0_rtv = 5.0e14_num
    REAL(num) :: tpfac, tpbar, tpmk, chi, alpha, ne
             
    tpfac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) tpfac = 2.0_num
    tpbar = tpfac*(gamma - 1.0_num)*(enbar - (1.0_num - nf)*ionise_pot)/(2.0_num - nf)             
    tpmk = tpbar*(Tp0/1.0e6_num)     

    ! a continuous radiation function requires slightly different
    ! chi values to those given in the LARE3D user guide
    IF (tpmk .GT. 0.02_num .AND. tpmk .LE. 0.0398_num) THEN
      chi = 1.2589_num
      alpha = 0.0_num
    ELSEIF (tpmk .LE. 0.0794_num) THEN
      chi = 794.328_num
      alpha = 2.0_num
    ELSEIF (tpmk .LE. 0.251_num) THEN
      chi = 5.01186_num
      alpha = 0.0_num
    ELSEIF (tpmk .LE. 0.562_num) THEN
      chi = 0.316228_num
      alpha = -2.0_num
    ELSEIF (tpmk .LE. 1.995_num) THEN
      chi = 1.0_num
      alpha = 0.0_num
    ELSEIF (tpmk .LE. 10.0_num) THEN
      chi = 1.58489_num
      alpha = -2.0_num/3.0_num
    ELSE
      chi = 0.0_num
      alpha = 0.0_num
    ENDIF

    ! electron number density in units of n0_rtv
    ne = (rhobar*(1.0_num - nf)*Dn0)/(m0*n0_rtv)
    
    ! return the radiative loss in SI units
    rad_loss = 2.87e-6_num*(ne**2)*chi*(tpmk**alpha)
    rad_alpha = alpha

  END SUBROUTINE calc_rtv_rad_loss


  ! calculate the radiation loss according to the method given by
  ! Klimchuk, J. A., Patsourakos, S., & Cargill, P. J., 2008, ApJ, 1351
  SUBROUTINE calc_kpc_rad_loss(rhobar, enbar, nf, rad_loss, rad_alpha)

    REAL(num), INTENT(IN) :: rhobar, enbar, nf
    REAL(num), INTENT(OUT) :: rad_loss, rad_alpha

    REAL(num) :: tpfac, tpbar, tp, tplog, chi, alpha, ne
    
    tpfac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) tpfac = 2.0_num
    tpbar = tpfac*(gamma - 1.0_num)*(enbar - (1.0_num - nf)*ionise_pot)/(2.0_num - nf)

    tp = tpbar*Tp0
    tplog = log10(tp)

    ! chi values have been converted from units of [erg s^-1 cm^3 K^-alpha] (see Klimchuk et al)
    ! to SI units of [J s^-1 m^3  K^-alpha], which simply requires multiplication by 1e-13
    IF (tp .LE. rad_min_tp) THEN
      chi = 0.0_num
      alpha = 0.0_num
    ELSEIF (tplog .LE. 4.97_num) THEN
      chi = 1.09e-44_num
      alpha = 2.0_num
    ELSEIF (tplog .LE. 5.67_num) THEN
      chi = 8.87e-30_num
      alpha = -1.0_num
    ELSEIF (tplog .LE. 6.18_num) THEN
      chi = 1.9e-35_num
      alpha = 0.0_num
    ELSEIF (tplog .LE. 6.55_num) THEN
      chi = 3.53e-26_num
      alpha = -1.5_num
    ELSEIF (tplog .LE. 6.9_num) THEN
      chi = 3.46e-38_num
      alpha = 1.0_num/3.0_num
    ELSEIF (tplog .LE. 7.63_num) THEN
      chi = 5.49e-29_num
      alpha = -1.0_num
    ELSE
      chi = 1.96e-40_num
      alpha = 0.5_num
    ENDIF

    ! electron number density in SI units
    ne = (rhobar*(1.0_num - nf)*Dn0)/m0

    ! return the radiative loss in SI units
    rad_loss = (ne**2)*chi*(tp**alpha)
    rad_alpha = alpha          

  END SUBROUTINE calc_kpc_rad_loss


  ! calculate the radiation loss according to the specified algorithm
  SUBROUTINE calc_rad_loss(rhobar, enbar, nf, rad_loss, rad_alpha)
    REAL(num), INTENT(IN) :: rhobar, enbar, nf
    REAL(num), INTENT(OUT) :: rad_loss, rad_alpha
    
    IF (rad_alg .EQ. 2) THEN
      CALL calc_kpc_rad_loss(rhobar, enbar, nf, rad_loss, rad_alpha)
    ELSE
      CALL calc_rtv_rad_loss(rhobar, enbar, nf, rad_loss, rad_alpha)
    ENDIF

    IF (rad_avTp) THEN
      rad_alpha = rad_alpha / 2.0_num
    ENDIF
  END SUBROUTINE calc_rad_loss
  

  ! return the timestep required to resolve radiation
  FUNCTION rad_dt    
    REAL(num) :: rad_dt
    REAL(num) :: rad_dt_local, loss_dt
    REAL(num) :: Lr, alpha

    rad_dt_local = largest_number  

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx

          CALL calc_rad_loss(rho(ix,iy,iz), energy(ix,iy,iz), xi_n(ix,iy,iz), Lr, alpha)
          
          IF (Lr .GE. rad_min_loss) THEN
            loss_dt = ABS((rad_dTp_frac*energy(ix,iy,iz)*En0*rho(ix,iy,iz)*Dn0) / (Lr*(1.0_num + rad_dTp_frac*alpha)))
            loss_dt = loss_dt/tm0
            rad_dt_local = MIN(rad_dt_local, loss_dt)
          ENDIF
  
        END DO
      END DO    		
    END DO

    rad_dt = rad_dt_local

  END FUNCTION rad_dt


  ! calculate the radiation losses explicitly
  SUBROUTINE radiate
    REAL(num) :: Lr, alpha
    REAL(num) :: Lr_enbar, a1, a2, rad_en
    REAL(num) :: rad_min_en

    rad_min_en = ((2.0_num/(gamma-1.0_num))*((kb_si*rad_min_tp)/m0))/En0

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx

          CALL calc_rad_loss(rho(ix,iy,iz), energy(ix,iy,iz), xi_n(ix,iy,iz), Lr, alpha)
            
          ! express Lr in terms of normalised energy units
          Lr_enbar = Lr*((dt*tm0)/(rho(ix,iy,iz)*Dn0))/En0
          a1 = alpha*Lr_enbar/energy(ix,iy,iz)
          a2 = Lr_enbar*(alpha - 1.0_num)

          rad_en = energy(ix,iy,iz) - (energy(ix,iy,iz) + a2)/(1.0_num + a1)
          IF (rad_en .GT. 0.0_num) THEN           
            IF (rad_min_en .GT. (energy(ix,iy,ix)-rad_en)) THEN
              rad_en = energy(ix,iy,iz) - rad_min_en
            END IF
            energy(ix,iy,iz) = energy(ix,iy,iz) - rad_en
            loss_rad = loss_rad + rad_en*rho(ix,iy,iz)*cv(ix,iy,iz)
          END IF
      
        END DO  
      END DO    		
    END DO

  END SUBROUTINE radiate


END MODULE conrad
