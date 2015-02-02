!*************************************************************************
! Controls all I/O and diagnostics. Output files are 'lare2d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'fort.5x'
! The idl package in 'plot.pro' gives simple loading and surface
! plotting based on these files. This isn't documented but is very simple!
!*************************************************************************
MODULE diagnostics

  USE shared_data
  USE boundary
  USE output_cartesian
  USE output
  USE iocontrol
  USE conrad
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines, energy_correction

CONTAINS

  SUBROUTINE output_routines(i) ! i = step index

    ! if halt set to false then code stops
    INTEGER, INTENT(IN) :: i

    INTEGER, PARAMETER :: out = 1000
    INTEGER, SAVE :: step = 1
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: data
    LOGICAL :: print_arrays, last_call
    REAL(num), DIMENSION(3) :: stagger = 0.0_num
    INTEGER, DIMENSION(3) :: dims
    
    ! this output routine uses the same sturcture as needed for mpi output
    ! this is more complicated than need for the serial code
    ! rank always equals zero in this serial code
    CHARACTER(LEN = 9+data_dir_max_length+n_zeros) :: filename
    CHARACTER(LEN = 35) :: filename_desc

    REAL(num) :: en_b, en_ke, en_int
    REAL(num) :: total_heating_visc, total_heating_ohmic, total_heating_dp
    REAL(num) :: total_htvisc_xy, total_htvisc_xz, total_htvisc_yz
    REAL(num) :: total_htvisc_xx, total_htvisc_yy, total_htvisc_zz
    REAL(num) :: total_con_supp, total_x_loss_con, total_y_loss_con, total_z_loss_con, total_loss_rad
    REAL(num) :: total_eta_crit_cnt
    REAL(num) :: max_jmag, max_temp
    REAL(num) :: total_rke_neg, total_rke_pos
            
    dims = (/ nx_global+1, ny_global+1, nz_global+1 /)
    
    IF (nsteps >= out) step = nsteps / out + 1 ! make sure output fits arrays
    
    IF (i == 0) THEN
      IF (restart) THEN
        CALL io_test(i, print_arrays, last_call)
        output_file = restart_snapshot+1
        RETURN
      ELSE
        IF (rank == 0) THEN
          CALL output_log
          WRITE(30) num, 23 
        END IF 
        output_file = 0
      END IF
    END IF

    IF (MOD(i, step) .EQ. 0 .OR. last_call) THEN ! do every (step) steps
      CALL energy_account(en_b, en_ke, en_int)
      
      CALL MPI_ALLREDUCE(heating_visc, total_heating_visc, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(heating_ohmic, total_heating_ohmic, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(heating_dp, total_heating_dp, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_xy, total_htvisc_xy, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_xz, total_htvisc_xz, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_yz, total_htvisc_yz, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_xx, total_htvisc_xx, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_yy, total_htvisc_yy, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(htvisc_zz, total_htvisc_zz, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(rke_neg, total_rke_neg, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(rke_pos, total_rke_pos, 1, mpireal, MPI_SUM, comm, errcode)
      
            
      IF (conduction) THEN
        IF (xbc_min .EQ. BC_OTHER .AND. (proc_x_min .EQ. MPI_PROC_NULL .OR. proc_x_max .EQ. MPI_PROC_NULL)) THEN
          x_loss_con = x_loss_con + con_x_loss()
        END IF
        IF (ybc_min .EQ. BC_OTHER .AND. (proc_y_min .EQ. MPI_PROC_NULL .OR. proc_y_max .EQ. MPI_PROC_NULL)) THEN
          y_loss_con = y_loss_con + con_y_loss()
        END IF
        IF (zbc_min .EQ. BC_OTHER .AND. (proc_z_min .EQ. MPI_PROC_NULL .OR. proc_z_max .EQ. MPI_PROC_NULL)) THEN
          z_loss_con = z_loss_con + con_z_loss()
        END IF
      END IF

      CALL MPI_ALLREDUCE(con_supp, total_con_supp, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(x_loss_con, total_x_loss_con, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(y_loss_con, total_y_loss_con, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(z_loss_con, total_z_loss_con, 1, mpireal, MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(loss_rad, total_loss_rad, 1, mpireal, MPI_SUM, comm, errcode)

      max_jmag = 0.0_num
      max_temp = 0.0_num                        
      CALL misc_account(max_jmag, max_temp)

      CALL MPI_ALLREDUCE(eta_crit_cnt, total_eta_crit_cnt, 1, mpireal, MPI_SUM, comm, errcode)      
      
      IF (rank .EQ. 0) THEN
        WRITE(30) time, en_b, en_ke, en_int
        WRITE(30) total_heating_visc, total_heating_ohmic, total_heating_dp
        WRITE(30) total_con_supp, total_x_loss_con, total_y_loss_con, total_z_loss_con, total_loss_rad
        WRITE(30) max_jmag, max_temp, total_eta_crit_cnt/(nx_global*ny_global*nz_global)
        WRITE(30) total_htvisc_xy, total_htvisc_xz, total_htvisc_yz, total_htvisc_xx, total_htvisc_yy, total_htvisc_zz
        WRITE(30) total_rke_neg, total_rke_pos        
      END IF      
    END IF

    CALL io_test(i, print_arrays, last_call) ! check if snapshot is needed
    
    IF (print_arrays) THEN ! output a snapshot of arrays
      IF (rank .EQ. 0) THEN
        WRITE(20, *) "Dumping ", output_file, " at time ", time, " (dt=", dt, ")" 
        CALL FLUSH(20)
      END IF

      ! Set the filename
      WRITE(filename_desc, '("(''nfs:'', a, ''/'', i", i3.3, ".", i3.3, &
          & ", ''.cfd'')")'), n_zeros, n_zeros
      WRITE(filename, filename_desc) TRIM(data_dir), output_file

      CALL cfd_open(filename, rank, comm, MPI_MODE_CREATE + MPI_MODE_WRONLY)
      CALL cfd_write_snapshot_data(REAL(time,dbl), i, 0)

      ALLOCATE(data(0:nx, 0:ny, 0:nz))
      
      CALL cfd_write_3d_cartesian_grid("Grid", "Grid", &
          xb_global(0:nx_global), yb_global(0:ny_global), &
          zb_global(0:nz_global), 0)
      
      IF (dump_mask(1))  THEN
        data = rho(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Rho", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(2))  THEN
        data = energy(0:nx, 0:ny, 0:nz)         
        CALL cfd_write_3d_cartesian_variable_parallel("Energy", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(3))  THEN
        data = vx(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Vx", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(4))  THEN
        data = vy(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Vy", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(5))  THEN
        data = vz(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Vz", "Velocity", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(6))  THEN
        data = bx(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Bx", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)         
      END IF

      IF (dump_mask(7))  THEN
        data = by(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("By", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(8))  THEN
        data = bz(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("Bz", "Magnetic_Field", &
            dims, stagger, "Grid", "Grid", data, subtype)        
      END IF

      IF (dump_mask(9)) THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
               data(ix,iy,iz) = (gamma - 1.0_num) &
                      * (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) &
                      / (2.0_num - xi_n(ix,iy,iz))                  
            END DO
          END DO
        END DO  
        
        CALL cfd_write_3d_cartesian_variable_parallel("Temperature", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(10))  THEN
        DO iz = 0, nz
          DO iy = 0, ny
            DO ix = 0, nx
               data(ix,iy,iz) = (energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz)) * ionise_pot) &
                          * (gamma - 1.0_num) * rho(ix,iy,iz)
            END DO
          END DO
        END DO
        CALL cfd_write_3d_cartesian_variable_parallel("Pressure", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(11)) THEN
        data = SQRT(gamma*(gamma - 1.0_num) * energy(0:nx,0:ny,0:nz)) 
        CALL cfd_write_3d_cartesian_variable_parallel("cs", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(12)) THEN
        data = parallel_current(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("j_par", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(13)) THEN
        data = perp_current(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("j_perp", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(14)) THEN
        data = xi_n(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("neutral_fraction", &
            "PIP", dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(15)) THEN
        data = eta_perp(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("eta_perp", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(16)) THEN
        data = eta(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("eta", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(17)) THEN
        data = jx_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("jx", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(18)) THEN
        data = jy_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("jy", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(19)) THEN
        data = jz_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("jz", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF      

      IF (dump_mask(20)) THEN
        data = p_visc(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("p_visc", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(21)) THEN
        data = visc_heat(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(22)) THEN
        data = visc_heat_xy(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_xy", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF      

      IF (dump_mask(23)) THEN
        data = visc_heat_xz(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_xz", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(24)) THEN
        data = visc_heat_yz(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_yz", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(25)) THEN
        data = visc_heat_xx(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_xx", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(26)) THEN
        data = visc_heat_yy(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_yy", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(27)) THEN
        data = visc_heat_zz(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("visc_heat_zz", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(28)) THEN
        data = ohmic_heat(0:nx, 0:ny, 0:nz)
        CALL cfd_write_3d_cartesian_variable_parallel("ohmic_heat", "Fluid", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF
      
      IF (dump_mask(29)) THEN
        data = fx_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("fx", "lorentz", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(30)) THEN
        data = fy_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("fy", "lorentz", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(31)) THEN
        data = fz_r(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("fz", "lorentz", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF
      
      IF (dump_mask(32)) THEN
        data = eta_num(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("eta_num", "PIP", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF

      IF (dump_mask(33)) THEN
        data = j_crit(0:nx, 0:ny, 0:nz) 
        CALL cfd_write_3d_cartesian_variable_parallel("j_crit", "current", &
            dims, stagger, "Grid", "Grid", data, subtype)
      END IF
      
      ! Close the file
      CALL cfd_close()

      DEALLOCATE(data)

      output_file = output_file + 1

    END IF           

    IF (last_call .AND. rank == 0) THEN ! output energy diagnostics etc
      WRITE(20, *) 'final nsteps / time =', i, time 
    END IF

  END SUBROUTINE output_routines



  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (restart) THEN
      ! avoid rewriting the snapshot you restarted from
      t1 = time + dt_snapshots
      restart = .FALSE.
    END IF

    IF (time >= t1) THEN
      print_arrays = .TRUE.
      t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test



  SUBROUTINE set_dt ! sets CFL limited step

    ! Assumes all variables are defined at the same point. Be careful
    ! with setting 'dt_multiplier' if you expect massive changes across
    ! cells.

    REAL(num) :: cons, dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt_local, dxlocal
    REAL(num) :: vxbp, vxbm, vybp, vybm, dvx, dvy, avxp, avxm, avyp, avym  
    REAL(num) :: vzbp, vzbm, dvz, avzp, avzm  
    REAL(num) :: dtr_local, cs, volume, ax, ay, az

    dt_local = largest_number
    dtr_local = largest_number
    cons = gamma * (gamma - 1.0_num)

    DO iz = -1, nz+2
      izm = iz - 1
      DO iy = -1, ny+2
        iym = iy - 1
        DO ix = -1, nx+2
          ixm = ix - 1

          w1 = bx(ix, iy, iz)**2 + by(ix, iy, iz)**2 + bz(ix, iy, iz)**2 
          cs = cons * energy(ix,iy,iz)      ! sound speed squared

          w2 = SQRT(cs + w1 / MAX(rho(ix, iy, iz), none_zero) &
              + 2.0_num * p_visc(ix, iy, iz) / MAX(rho(ix, iy, iz), none_zero)) 

          ! find ideal MHD CFL limit for Lagrangian step
          dt1 = MIN(dxb(ix), dyb(iy), dzb(iz)) / w2 
          dt_local = MIN(dt_local, dt1)

          ax = 0.25_num * dyb(iy) * dzb(iz) 
          vxbp = (vx(ix,iy,iz) + vx(ix,iym,iz) + vx(ix,iy,izm) + vx(ix,iym,izm)) * ax
          vxbm = (vx(ixm,iy,iz) + vx(ixm,iym,iz) + vx(ixm,iy,izm) + vx(ixm,iym,izm)) * ax
          ay = 0.25_num * dxb(ix) * dzb(iz) 
          vybp = (vy(ix,iy,iz) + vy(ixm,iy,iz) + vy(ix,iy,izm) + vy(ixm,iy,izm)) * ay
          vybm = (vy(ix,iym,iz) + vy(ixm,iym,iz) + vy(ix,iym,izm) + vy(ixm,iym,izm)) * ay
          az = 0.25_num * dyb(iy) * dxb(ix) 
          vzbp = (vz(ix,iy,iz) + vz(ix,iym,iz) + vz(ixm,iy,iz) + vz(ixm,iym,iz)) * az
          vzbm = (vz(ix,iy,izm) + vz(ix,iym,izm) + vz(ixm,iy,izm) + vz(ixm,iym,izm)) * az
          
          dvx = ABS(vxbp - vxbm)
          dvy = ABS(vybp - vybm) 
          dvz = ABS(vzbp - vzbm) 
          avxp = ABS(vxbp)
          avxm = ABS(vxbm)
          avyp = ABS(vybp)
          avym = ABS(vybm)
          avzp = ABS(vzbp)
          avzm = ABS(vzbm)
           
          volume = ax * dxb(ix)
          dt5 = volume / MAX(avxp, avxm, dvx, 1.e-10_num * volume)
          dt6 = volume / MAX(avyp, avym, dvy, 1.e-10_num  * volume)
          dt7 = volume / MAX(avzp, avzm, dvz, 1.e-10_num  * volume)
          
          ! fix dt for remap step 
          dt_local = MIN(dt_local, dt5, dt6, dt7)

          ! note resistive limits assumes uniform resistivity hence cautious
          ! factor 0.2
          dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 &
              + 1.0_num / dyb(iy)**2 + 1.0_num / dzb(iz)**2)

          IF (cowling_resistivity) THEN
            dt3 = 0.2_num * dxlocal &
                / MAX(MAX(eta(ix, iy, iz), eta_perp(ix, iy, iz)), none_zero)
          ELSE
            dt3 = 0.2_num * dxlocal / MAX(eta(ix, iy, iz), none_zero)
          END IF

          ! adjust to accomodate resistive effects
          dtr_local = MIN(dtr_local, dt3)

        END DO
      END DO
    END DO

    IF (radiation .AND. (.NOT. rad_implicit)) THEN
      dt_local = MIN(dt_local, rad_dt())
    END IF    

    CALL MPI_ALLREDUCE(dt_local, dt, 1, mpireal, MPI_MIN, comm, errcode)
    CALL MPI_ALLREDUCE(dtr_local, dtr, 1, mpireal, MPI_MIN, comm, errcode)

    dtr = dt_multiplier * dtr
    dt = dt_multiplier * dt

    time = time + dt

    IF (rank == 0 .AND. radiation .AND. (.NOT. rad_implicit)) THEN
      PRINT * , time, ": dt = ", dt, "."
    END IF

  END SUBROUTINE set_dt



  SUBROUTINE energy_account(energy_b, energy_ke, energy_int)

    REAL(num), INTENT(OUT) :: energy_b, energy_ke, energy_int
    REAL(dbl) :: energy_b_local, energy_ke_local, energy_int_local
    REAL(num) :: cv_v, rho_v, a, b, c
    
    energy_b_local = 0.0_dbl
    energy_ke_local = 0.0_dbl
    energy_int_local = 0.0_dbl

    DO iz = 1, nz
      izm = iz - 1
      DO iy = 1, ny
        iym = iy - 1
        DO ix = 1, nx
          ixm = ix - 1

          w2 = (bx(ix, iy, iz)**2 + bx(ixm, iy, iz)**2) / 2.0_dbl
          w3 = (by(ix, iy, iz)**2 + by(ix, iym, iz)**2) / 2.0_dbl
          w4 = (bz(ix, iy, iz)**2 + bz(ix, iy, izm)**2) / 2.0_dbl
          w1 = (w2 + w3 + w4) / 2.0_dbl
          energy_b_local = energy_b_local + w1 * cv(ix, iy, iz)

          energy_int_local = energy_int_local &
              + energy(ix, iy, iz) * rho(ix, iy, iz) * cv(ix, iy, iz)
                                  
        END DO
      END DO
    END DO

    DO iz = 0, nz
      izp = iz + 1
      DO iy = 0, ny
        iyp = iy + 1
        DO ix = 0, nx
          ixp = ix + 1

          ! WARNING the KE is summed on the vertices
          rho_v = (rho(ix, iy, iz) * cv(ix, iy, iz) &
                + rho(ixp, iy , iz ) * cv(ixp, iy , iz ) &
                + rho(ix , iyp, iz ) * cv(ix , iyp, iz ) &
                + rho(ixp, iyp, iz ) * cv(ixp, iyp, iz ) &
                + rho(ix , iy , izp) * cv(ix , iy , izp) &
                + rho(ixp, iy , izp) * cv(ixp, iy , izp) &
                + rho(ix , iyp, izp) * cv(ix , iyp, izp) &
                + rho(ixp, iyp, izp) * cv(ixp, iyp, izp))

          cv_v = (cv(ix, iy, iz) + cv(ixp, iy, iz) &
               + cv(ix, iyp, iz ) + cv(ixp, iyp, iz ) &
               + cv(ix, iy , izp) + cv(ixp, iy , izp) &
               + cv(ix, iyp, izp) + cv(ixp, iyp, izp))

          rho_v = rho_v / cv_v
          cv_v = cv_v / 8.0_dbl
          w1 = rho_v * cv_v * (vx(ix, iy, iz)**2 + vy(ix, iy, iz)**2 + vz(ix, iy, iz)**2)

          IF ((ix == 0) .OR. (ix == nx)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          IF ((iy == 0) .OR. (iy == ny)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          IF ((iz == 0) .OR. (iz == nz)) THEN
            w1 = w1 / 2.0_dbl
          END IF

          energy_ke_local = energy_ke_local + w1 / 2.0_dbl

        END DO
      END DO
    END DO       

    a = REAL(energy_ke_local, num)
    b = REAL(energy_b_local, num)
    c = REAL(energy_int_local, num)
    
    CALL MPI_ALLREDUCE(a, energy_ke, 1, mpireal, MPI_SUM, comm, errcode)
    CALL MPI_ALLREDUCE(b, energy_b, 1, mpireal, MPI_SUM, comm, errcode)
    CALL MPI_ALLREDUCE(c, energy_int, 1, mpireal, MPI_SUM, comm, errcode)
    
  END SUBROUTINE energy_account



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    !WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    !delta_ke = delta_ke / (rho * cv)

    rke_neg = 0.0_num
    rke_pos = 0.0_num
        
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
            IF (delta_ke(ix,iy,iz) .GT. 0.0_num) THEN
              rke_neg = rke_neg + delta_ke(ix,iy,iz)
              energy(ix,iy,iz) = energy(ix,iy,iz) + delta_ke(ix,iy,iz)/(rho(ix,iy,iz)*cv(ix,iy,iz))              
            ELSE
              rke_pos = rke_pos + ABS(delta_ke(ix,iy,iz))
              delta_ke(ix,iy,iz) = 0.0_num
            END IF    
        END DO
      END DO
    END DO            

    CALL energy_bcs

  END SUBROUTINE energy_correction


  SUBROUTINE misc_account(max_jmag, max_temp)

    REAL(num), INTENT(OUT) :: max_jmag, max_temp
                
    REAL(num) :: jxe, jye, jze, jxpe, jype, jzpe, jxv, jyv, jzv
    REAL(num) :: temp_fac, jmag, temp

    temp_fac = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. neutral_gas) temp_fac = 2.0_num

    
    ! calculate jmag and temp
    jmag = 0.0_num 
    temp = 0.0_num
  
    DO ix = 0, nx
      DO iy = 0, ny
        DO iz = 0, nz
        
          ixp = ix + 1
          iyp = iy + 1
          izp = iz + 1  
                 
          ! edge-centred currents 
          jxe = (bz(ix,iyp,iz) - bz(ix,iy,iz))/dyc(iy) - (by(ix,iy,izp) - by(ix,iy,iz))/dzc(iz)       
          jye = (bx(ix,iy,izp) - bx(ix,iy,iz))/dzc(iz) - (bz(ixp,iy,iz) - bz(ix,iy,iz))/dxc(ix)
          jze = (by(ixp,iy,iz) - by(ix,iy,iz))/dxc(ix) - (bx(ix,iyp,iz) - bx(ix,iy,iz))/dyc(iy)

          ! vertex-centred currents next cell on
          jxpe = (bz(ixp,iyp,iz) - bz(ixp,iy,iz))/dyc(iy) - (by(ixp,iy,izp) - by(ixp,iy,iz))/dzc(iz)               
          jype = (bx(ix,iyp,izp) - bx(ix,iyp,iz))/dzc(iz) - (bz(ixp,iyp,iz) - bz(ix,iyp,iz))/dxc(ix)               
          jzpe = (by(ixp,iy,izp) - by(ix,iy,izp))/dxc(ix) - (bx(ix,iyp,izp) - bx(ix,iy,izp))/dyc(iy)

          ! vertex-centred currents
          jxv = 0.5_num * (jxe + jxpe) 
          jyv = 0.5_num * (jye + jype) 
          jzv = 0.5_num * (jze + jzpe)

          jmag = MAX(SQRT(jxv**2 + jyv**2 + jzv**2), jmag)
          temp = MAX(temp_fac*(gamma - 1.0_num)*(energy(ix,iy,iz) - (1.0_num - xi_n(ix,iy,iz))*ionise_pot)/(2.0_num - xi_n(ix,iy,iz)), temp)                  
                          
        ENDDO        
      ENDDO
    ENDDO
    
    CALL MPI_ALLREDUCE(jmag, max_jmag, 1, mpireal, MPI_MAX, comm, errcode)
    CALL MPI_ALLREDUCE(temp, max_temp, 1, mpireal, MPI_MAX, comm, errcode)      

  END SUBROUTINE misc_account


  SUBROUTINE output_log ! writes basic data to 'lare3d.dat'

    WRITE(20, *) 'nprocx, nprocy, nproca = ', nprocx, nprocy, nprocz
    WRITE(20, *) 'nx, ny, nz = ', nx, ny, nz
    WRITE(20, *)
    WRITE(20, *) 'length_x = ', length_x 
    WRITE(20, *) 'length_y = ', length_y 
    WRITE(20, *) 'length_z = ', length_z
    WRITE(20, *)
    WRITE(20, *) 'x_start = ', x_start
    WRITE(20, *) 'x_end = ', x_end
    WRITE(20, *) 'x_stretch = ', x_stretch
    WRITE(20, *) 'y_start = ', y_start
    WRITE(20, *) 'y_end = ', y_end
    WRITE(20, *) 'y_stretch = ', y_stretch
    WRITE(20, *) 'z_start = ', z_start
    WRITE(20, *) 'z_end = ', z_end
    WRITE(20, *) 'z_stretch = ', z_stretch
    WRITE(20, *)
    WRITE(20, *) 'xbc_min = ', xbc_min
    WRITE(20, *) 'xbc_max = ', xbc_max
    WRITE(20, *) 'ybc_max = ', ybc_min
    WRITE(20, *) 'ybc_min = ', ybc_max
    WRITE(20, *) 'zbc_min = ', zbc_min
    WRITE(20, *) 'zbc_max = ', zbc_max
    WRITE(20, *)
    WRITE(20, *) 't_start, t_end = ', time, t_end
    WRITE(20, *) 'nsteps =', nsteps
    WRITE(20, *) 'dt_multiplier =', dt_multiplier
    WRITE(20, *)
    WRITE(20, *)
    WRITE(20, *) 'resistive_mhd = ', resistive_mhd   
    WRITE(20, *) 'eta_background = ', eta_background
    WRITE(20, *) 'eta_crit = ', eta_crit
    WRITE(20, *)
    WRITE(20, *) 'conduction = ', conduction 
    WRITE(20, *) 'con_flux_limited = ', con_flux_limited
    WRITE(20, *) 'con_limit = ', con_limit
    WRITE(20, *) 'con_b_min = ', con_b_min
    WRITE(20, *)
    WRITE(20, *) 'radiation = ', radiation
    WRITE(20, *) 'rad_implicit = ', rad_implicit
    WRITE(20, *) 'rad_alg = ', rad_alg 
    WRITE(20, *) 'rad_min_tp = ', rad_min_tp
    WRITE(20, *) 'rad_min_loss = ', rad_min_loss 
    WRITE(20, *) 'rad_dTp_frac = ', rad_dTp_frac
    WRITE(20, *) 'rad_avTp = ', rad_avTp
    WRITE(20, *)
    WRITE(20, *) 'sor_it_lim = ', sor_it_lim    
    WRITE(20, *) 'sor_w = ', sor_w    
    WRITE(20, *) 'sor_err = ', sor_err
    WRITE(20, *) 'sor_err_wgt = ', sor_err_wgt
    WRITE(20, *)
    WRITE(20, *) 'conrad_half_dt = ', conrad_half_dt
    WRITE(20, *)
    WRITE(20, *) 'rke = ', rke
    WRITE(20, *)
    WRITE(20, *) 'initial = ', initial
    WRITE(20, *) 'restart_snapshot = ', restart_snapshot
    WRITE(20, *) 'cowling_resistivity = ', cowling_resistivity
    WRITE(20, *) 'ionise_pot = ', ionise_pot
    WRITE(20, *)
#ifndef Q_MONO
    WRITE(20, *) 'tensor shock viscosity'
#else
    WRITE(20, *) 'q_mono viscosity'
#endif
    WRITE(20, *) 'linear viscosity coeff = ', visc1
    WRITE(20, *) 'quadratic viscosity coeff = ', visc2
    WRITE(20, *) 'uniform tensor viscosity coeff = ', visc3 
    WRITE(20, *)
    WRITE(20, *) 'damping = ', damping
    WRITE(20, *) 'eos_number = ', eos_number
    WRITE(20, *) 'neutral_gas = ', neutral_gas
    WRITE(20, *)
    WRITE(20, *) 'iniTp = ', iniTp
    WRITE(20, *) 'iniEn = ', iniEn
    WRITE(20, *)
    WRITE(20, *) 'lp_rad = ', lp_rad
    WRITE(20, *) 'lp_alp = ', lp_alp
    WRITE(20, *) 'lp_omg = ', lp_omg
    WRITE(20, *)
    WRITE(20, *) 'perturbation = ', perturbation
    WRITE(20, *) 'ptb_sign = ', ptb_sign
    WRITE(20, *) 'ptb_k = ', ptb_k
    WRITE(20, *) 'ptb_amp = ', ptb_amp
    WRITE(20, *) 'ptb_lim = ', ptb_lim
    WRITE(20, *)
    WRITE(20, *) 'diag_corks = ', diag_corks
    WRITE(20, *)
    WRITE(20, *) 'gamma = ', gamma
    WRITE(20, *) 'mf = ', mf
    WRITE(20, *) 'm0 = ', m0
    WRITE(20, *)
    WRITE(20, *) 'B0 = ', B0
    WRITE(20, *) 'L0 = ', L0
    WRITE(20, *) 'NmDn0 = ', NmDn0
    WRITE(20, *)
    WRITE(20, *) 'Dn0  = ', Dn0 
    WRITE(20, *) 'vA = ', vA
    WRITE(20, *) 'tm0 = ', tm0
    WRITE(20, *) 'P0 = ', P0
    WRITE(20, *) 'Tp0 = ', Tp0
    WRITE(20, *) 'Tpbar = ', iniTp/Tp0
    WRITE(20, *) 'En0 = ', En0
    WRITE(20, *) 'Enbar = ', iniEn/En0
    WRITE(20, *) 'K0 = ', K0
    WRITE(20, *) 'eta0 = ', eta0
    
  END SUBROUTINE output_log

END MODULE diagnostics
