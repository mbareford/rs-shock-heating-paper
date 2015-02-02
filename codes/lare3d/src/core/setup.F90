MODULE setup

  USE shared_data
  USE iocommon
  USE iocontrol
  USE input
  USE input_cartesian

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: before_control, after_control
  PUBLIC :: grid
  PUBLIC :: open_files, close_files, restart_data

  REAL(num), DIMENSION(:), ALLOCATABLE :: dxnew, dynew, dznew

CONTAINS

  SUBROUTINE before_control
    ! Setup basic variables which have to have default values

    nprocx = 0
    nprocy = 0
    nprocz = 0

    time = 0.0_num
    gamma = 5.0_num / 3.0_num

    IF (num .EQ. 4) mpireal = MPI_REAL

  END SUBROUTINE before_control



  SUBROUTINE after_control
    ! Setup arrays and other variables which can only be set after
    ! user input

    p_visc = 0.0_num
    visc_heat = 0.0_num
    ohmic_heat = 0.0_num

    visc_heat_xy = 0.0_num
    visc_heat_xz = 0.0_num
    visc_heat_yz = 0.0_num
    visc_heat_xx = 0.0_num
    visc_heat_yy = 0.0_num
    visc_heat_zz = 0.0_num

    dis_dzbx_w = 0.0_num
    dis_dzby_w = 0.0_num

    eta = 0.0_num
    eta_num = 0.0_num
    j_crit = 0.0_num
    grav = 0.0_num

    rho = 0.0_num
    energy = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

  END SUBROUTINE after_control



  SUBROUTINE grid ! stretched and staggered grid

    REAL(num) :: dx, dy, dz, xcstar, ycstar, zcstar
    INTEGER :: ix, iy, n0, n1
    INTEGER :: nx0, ny0, nz0
    INTEGER :: nxp, nyp, nzp
    INTEGER :: cx, cy, cz

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx .NE. nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
    ELSE
      nxp = nprocx
    ENDIF
    IF (ny0 * nprocy .NE. ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
    ELSE
      nyp = nprocy
    ENDIF
    IF (nz0 * nprocz .NE. nz_global) THEN
      nzp = (nz0 + 1) * nprocz - nz_global
    ELSE
      nzp = nprocz
    ENDIF

    ALLOCATE(dxnew(-2:nx_global+2))
    ALLOCATE(dynew(-2:ny_global+2))
    ALLOCATE(dznew(-2:nz_global+2))

    ! initially assume uniform grid
    dx = 1.0_num / REAL(nx_global, num)
    dy = 1.0_num / REAL(ny_global, num)
    dz = 1.0_num / REAL(nz_global, num)

    length_x = x_end - x_start
    length_y = y_end - y_start
    length_z = z_end - z_start

    xb_global(0) = 0.0_num ! grid cell boundary for x coordinates
    DO ix = -2, nx_global+2
      xb_global(ix) = xb_global(0) + REAL(ix, num) * dx
    END DO
    xb_global = xb_global * (x_end - x_start) + x_start

    IF (x_stretch) CALL stretch_x ! stretch grid ?

    ! define position of ghost cells using sizes of adjacent cells
    xb_global(nx_global+1) = xb_global(nx_global) &
        + (xb_global(nx_global) - xb_global(nx_global-1))

    ! needed for ghost cell
    xb_global(nx_global+2) = xb_global(nx_global+1) &
        + (xb_global(nx_global+1) - xb_global(nx_global))
    xb_global(-1) = xb_global( 0) - (xb_global(1) - xb_global( 0))
    xb_global(-2) = xb_global(-1) - (xb_global(0) - xb_global(-1))

    cx = coordinates(3)
    IF (cx .LT. nxp) THEN
      n0 = cx * nx0
      n1 = (cx + 1) * nx0
    ELSE
      n0 = nxp * nx0 + (cx - nxp) * (nx0 + 1)
      n1 = nxp * nx0 + (cx - nxp + 1) * (nx0 + 1)
    ENDIF

    xb = xb_global(n0-2:n1+2)

    DO ix = -1, nx+2
      ixm = ix - 1
      xc(ix) = 0.5_num * (xb(ixm) + xb(ix)) ! cell centre
    END DO

    DO ix = -1, nx+1
      ixp = ix + 1
      dxc(ix) = xc(ixp) - xc(ix) ! distance between centres
    END DO

    IF (cx == nprocx - 1) THEN
      dxc(nx+2) = dxc(nx+1)
    ELSE
      xcstar = 0.5_num * (xb(nx+2) + xb_global(n1+3))
      dxc(nx+2) = xcstar - xc(nx+2)
    END IF

    DO ix = -1, nx+2
      ixm = ix - 1
      dxb(ix) = xb(ix) - xb(ixm) ! cell width
    END DO

    yb_global(0) = 0.0_num ! grid cell boundary for y coordinates
    DO iy = -2, ny_global+2
      yb_global(iy) = yb_global(0) + REAL(iy, num) * dy
    END DO
    yb_global = yb_global * (y_end - y_start) + y_start

    IF (y_stretch) CALL stretch_y ! stretch grid ?

    ! define position of ghost cells using sizes of adjacent cells
    yb_global(ny_global+1) = yb_global(ny_global) &
        + (yb_global(ny_global) - yb_global(ny_global-1))

    ! needed for ghost cell
    yb_global(ny_global+2) = yb_global(ny_global+1) &
        + (yb_global(ny_global+1) - yb_global(ny_global))
    yb_global(-1) = yb_global( 0) - (yb_global(1) - yb_global( 0))
    yb_global(-2) = yb_global(-1) - (yb_global(0) - yb_global(-1))

    cy = coordinates(2)
    IF (cy .LT. nyp) THEN
      n0 = cy * ny0
      n1 = (cy + 1) * ny0
    ELSE
      n0 = nyp * ny0 + (cy - nyp) * (ny0 + 1)
      n1 = nyp * ny0 + (cy - nyp + 1) * (ny0 + 1)
    ENDIF

    yb = yb_global(n0-2:n1+2)

    DO iy = -1, ny+2
      iym = iy - 1
      yc(iy) = 0.5_num * (yb(iym) + yb(iy)) ! cell centre
    END DO

    DO iy = -1, ny+1
      iyp = iy + 1
      dyc(iy) = yc(iyp) - yc(iy) ! distance between centres
    END DO

    IF (cy == nprocy - 1) THEN
      dyc(ny+2) = dyc(ny+1)
    ELSE
      ycstar = 0.5_num * (yb(ny+2) + yb_global(n1+3))
      dyc(ny+2) = ycstar - yc(ny+2)
    END IF

    DO iy = -1, ny+2
      iym = iy - 1
      dyb(iy) = yb(iy) - yb(iym) ! cell width
    END DO

    zb_global(0) = 0.0_num ! grid cell boundary for z coordinates
    DO iz = -2, nz_global+2
      zb_global(iz) = zb_global(0) + REAL(iz, num) * dz
    END DO
    zb_global = zb_global * (z_end - z_start) + z_start

    IF (z_stretch) CALL stretch_z ! stretch grid ?

    ! define position of ghost cells using sizes of adjacent cells
    zb_global(nz_global+1) = zb_global(nz_global) &
        + (zb_global(nz_global) - zb_global(nz_global-1))

    ! needed for ghost cell
    zb_global(nz_global+2) = zb_global(nz_global+1) &
        + (zb_global(nz_global+1) - zb_global(nz_global))

    zb_global(-1) = zb_global( 0) - (zb_global(1) - zb_global( 0))
    zb_global(-2) = zb_global(-1) - (zb_global(0) - zb_global(-1))

    cz = coordinates(1)
    IF (cz .LT. nzp) THEN
      n0 = cz * nz0
      n1 = (cz + 1) * nz0
    ELSE
      n0 = nzp * nz0 + (cz - nzp) * (nz0 + 1)
      n1 = nzp * nz0 + (cz - nzp + 1) * (nz0 + 1)
    ENDIF

    zb = zb_global(n0-2:n1+2)

    DO iz = -1, nz+2
      izm = iz - 1
      zc(iz) = 0.5_num * (zb(izm) + zb(iz)) ! cell centre
    END DO

    DO iz = -1, nz+1
      izp = iz + 1
      dzc(iz) = zc(izp) - zc(iz) ! distance between centres
    END DO

    IF (cz == nprocz - 1) THEN
      dzc(nz+2) = dzc(nz+1)
    ELSE
      zcstar = 0.5_num * (zb(nz+2) + zb_global(n1+3))
      dzc(nz+2) = zcstar - zc(nz+2)
    END IF

    DO iz = -1, nz+2
      izm = iz - 1
      dzb(iz) = zb(iz) - zb(izm) ! cell width
    END DO

    DO ix = -1, nx+2
      DO iy = -1, ny+2
        DO iz = -1, nz+2
          cv(ix, iy, iz) = dxb(ix) * dyb(iy) * dzb(iz) ! define the cell area
        END DO
      END DO
    END DO

    DEALLOCATE(dxnew, dynew, dznew)

  END SUBROUTINE grid



  ! Subroutine stretches the grid in the x direction
  SUBROUTINE stretch_x ! replace with any stretching algorithm as needed

    REAL(num) :: width, dx, L, f, lx_new

    ! new total length
    lx_new = 200.0_num

    ! centre of tanh stretching in unstretched coordinates
    L = length_x / 1.5_num

    ! width of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num

    f = (lx_new - length_x) / (length_x - L) / 2.0_num

    dx = length_x / REAL(nx_global, num)
    dxnew = dx + f * (1.0_num + TANH((ABS(xb_global) - L) / width)) * dx

!!$    DO ix = nx_global/2+1, nx_global+2
!!$      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
!!$    END DO
!!$    DO ix = nx_global/2-1, -2, -1
!!$      xb_global(ix) = xb_global(ix+1) - dxnew(ix)
!!$    END DO

    DO ix = 1, nx_global+2
      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    END DO

    length_x = lx_new

  END SUBROUTINE stretch_x



  ! Subroutine stretches the domain in the y direction
  SUBROUTINE stretch_y ! stretch domain upwards only

    REAL(num) :: width, dy, L, f, ly_new

    ! new tolal length
    ly_new = 100.0_num

    ! centre of tanh stretching in unstretched coordinates
    L = length_y / 1.5_num

    ! width of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num

    f = (ly_new - length_y) / (length_y - L) / 2.0_num

    dy = length_y / REAL(ny_global, num)
    dynew = dy + f * (1.0_num + TANH((ABS(yb_global) - L) / width)) * dy

!!$    DO iy = ny_global/2+1, ny_global+2
!!$      yb_global(iy) = yb_global(iy-1) + dynew(iy)
!!$    END DO
!!$    DO iy = ny_global/2-1, -2, -1
!!$      yb_global(iy) = yb_global(iy+1) - dynew(iy)
!!$    END DO

    DO iy = 1, ny_global+2
      yb_global(iy) = yb_global(iy-1) + dynew(iy)
    END DO

    length_y = ly_new

  END SUBROUTINE stretch_y



  ! Subroutine stretches the domain in the z direction
  SUBROUTINE stretch_z

    REAL(num) :: width, dz, L, f, lz_new

    ! new tolal length
    lz_new = 33.0_num

    ! centre of tanh stretching in unstretched coordinates
    L = 2.0_num * length_z / 3.0_num

    ! width of tanh stretching in unstretched coordinates
    width = length_z / 10.0_num

    f = (lz_new - length_z) / (length_z - L) / 2.0_num

    dz = length_z / REAL(nz_global, num)
    dznew = dz + f * (1.0_num + TANH((ABS(zb_global) - L) / width)) * dz

    DO iz = 1, nz_global+2
      zb_global(iz) = zb_global(iz-1) + dznew(iz)
    END DO

  END SUBROUTINE stretch_z



  ! Open the output diagnostic files
  SUBROUTINE open_files

    CHARACTER(LEN = 11+data_dir_max_length) :: file2
    CHARACTER(LEN = 7+data_dir_max_length) :: file3
    CHARACTER(LEN = 8+data_dir_max_length) :: file4
    INTEGER :: ios, dum_int1, dum_int2
    REAL(num) :: tm, dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8

    IF (rank == 0) THEN
      WRITE(file2, '(a, "/lare3d.dat")') TRIM(data_dir)
      IF (IAND(initial, IC_RESTART) .NE. 0) THEN
        OPEN(unit = 20, STATUS = 'OLD', FILE = file2, POSITION = 'APPEND', iostat = ios)
      ELSE
        OPEN(unit = 20, STATUS = 'NEW', FILE = file2, iostat = ios)
      END IF

      IF (ios .NE. 0) THEN
        PRINT *, "Unable to open file lare3d.dat for writing. This is ", &
                 "most commonly caused by the output directory not existing"
        PRINT *, " "
        PRINT *, " "
        CALL MPI_ABORT(comm, errcode)
      END IF

      
      IF (IAND(initial, IC_RESTART) .NE. 0) THEN
        WRITE(file3, '(a, "/en.dat")') TRIM(data_dir)
        WRITE(file4, '(a, "/en1.dat")') TRIM(data_dir)    
        OPEN(unit=40, STATUS='OLD', FILE=file3, FORM='UNFORMATTED', ACCESS='STREAM', ACTION='READ', iostat=ios)
        OPEN(unit=30, STATUS='NEW', FILE=file4, FORM='UNFORMATTED', ACCESS='STREAM', iostat=ios)
        READ(unit=40, iostat=ios, ERR=901, END=902) dum_int1, dum_int2
        WRITE(30) dum_int1, dum_int2
        DO 
          READ(unit=40, iostat=ios, ERR=901, END=902) tm, dum1, dum2, dum3, &
                                                      heating_visc, heating_ohmic, heating_dp, &
                                                      con_supp, x_loss_con, y_loss_con, z_loss_con, loss_rad, &
                                                      dum4, dum5, dum6, &
                                                      htvisc_xy, htvisc_xz, htvisc_yz, htvisc_xx, htvisc_yy, htvisc_zz, &
                                                      dum7, dum8 
        
          IF (tm .LE. time) THEN
            WRITE(30) tm, dum1, dum2, dum3
            WRITE(30) heating_visc, heating_ohmic, heating_dp
            WRITE(30) con_supp, x_loss_con, y_loss_con, z_loss_con, loss_rad
            WRITE(30) dum4, dum5, dum6
            WRITE(30) htvisc_xy, htvisc_xz, htvisc_yz, htvisc_xx, htvisc_yy, htvisc_zz
            WRITE(30) dum7, dum8
          ELSE
            GOTO 902
          ENDIF
               
        WRITE(30) 
        WRITE(30) 

        ENDDO
901     PRINT *, "Error ", ios, " when reading en.dat file in preparation for restart."
902     CLOSE(unit=40)

        PRINT *, "Restart value for heating_visc is ", heating_visc, "."
        PRINT *, "Restart value for heating_ohmic is ", heating_ohmic, "."
        PRINT *, "Restart value for heating_dp is ", heating_dp, "."
        PRINT *, "Restart value for con_supp is ", con_supp, "."
        PRINT *, "Restart value for x_loss_con is ", x_loss_con, "."
        PRINT *, "Restart value for y_loss_con is ", y_loss_con, "."
        PRINT *, "Restart value for z_loss_con is ", z_loss_con, "."
        PRINT *, "Restart value for loss_rad is ", loss_rad, "."  
        PRINT *, "Restart value for htvisc_xy is ", htvisc_xy, "."  
        PRINT *, "Restart value for htvisc_xz is ", htvisc_xz, "."
        PRINT *, "Restart value for htvisc_yz is ", htvisc_yz, "."
        PRINT *, "Restart value for htvisc_xx is ", htvisc_xx, "."
        PRINT *, "Restart value for htvisc_yy is ", htvisc_yy, "."  
        PRINT *, "Restart value for htvisc_zz is ", htvisc_zz, "."          

      ELSE
        WRITE(file3, '(a, "/en.dat")') TRIM(data_dir)
        OPEN(unit=30, STATUS='NEW', FILE=file3, FORM='UNFORMATTED', ACCESS='STREAM', iostat=ios)
      END IF

      IF (ios .NE. 0) THEN
        PRINT *, "Unable to open file en.dat for writing. This is ", &
                 "most commonly caused by the output directory not existing"
        PRINT *, " "
        PRINT *, " "
        CALL MPI_ABORT(comm, errcode)
      END IF
    END IF

  END SUBROUTINE open_files



  ! Close the output diagnostic files
  SUBROUTINE close_files

    IF (rank == 0) THEN
      CLOSE(unit = 20)
      CLOSE(unit = 30)
    END IF

  END SUBROUTINE close_files



  ! Subroutine to perform string comparisons
  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) :: str_in, str_test
    CHARACTER(30) :: str_trim
    LOGICAL :: str_cmp

    str_trim = TRIM(ADJUSTL(str_in))

    IF (LEN(str_test) .GT. LEN(str_in)) THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    IF (str_trim(LEN(str_test)+1:LEN(str_test)+1) .NE. " ") THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    str_cmp = str_trim(1:LEN(str_test)) == str_test

  END FUNCTION str_cmp



  ! Restart from previous output dumps
  SUBROUTINE restart_data

    CHARACTER(LEN = 20+data_dir_max_length) :: filename
    CHARACTER(LEN = 20) :: name, class, mesh_name, mesh_class
    INTEGER :: nblocks, type, nd, sof, snap
    INTEGER, DIMENSION(3) :: dims
    REAL(dbl) :: time_d
    REAL(num), DIMENSION(3) :: extent
    REAL(num), DIMENSION(3) :: stagger
    REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: data
    REAL(num) :: ix2,iy2,iz2

    ! Create the filename for the last snapshot
#ifdef MHDCLUSTER
    WRITE(filename, '("nfs:", a, "/", i4.4, ".cfd")') &
        TRIM(data_dir), restart_snapshot
#else
    WRITE(filename, '(a, "/", i4.4, ".cfd")') &
        TRIM(data_dir), restart_snapshot
#endif

    output_file = restart_snapshot


    ALLOCATE(data(0:nx, 0:ny, 0:nz))
    CALL cfd_open(filename, rank, comm, MPI_MODE_RDONLY)
    ! Open the file
    nblocks = cfd_get_nblocks()

    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, type)
      IF (rank == 0) PRINT *, ix, name, class, type

      IF (type == TYPE_SNAPSHOT) THEN 
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
      END IF

      IF (type == TYPE_MESH) THEN
        ! Strangely, LARE doesn't actually read in the grid from a file
        ! This can be fixed, but for the moment, just go with the flow and
        ! Replicate the old behaviour

        CALL cfd_skip_block()
      ELSE IF (type == TYPE_MESH_VARIABLE) THEN
        CALL cfd_get_common_meshtype_metadata_all(type, nd, sof)

        IF (nd /= DIMENSION_3D) THEN
          IF (rank == 0) PRINT *, "Non 3D Dataset found in input file, ", &
              "ignoring and continuting."
          CALL cfd_skip_block()
          CYCLE
        END IF

        IF (type /= VAR_CARTESIAN) THEN
          IF (rank == 0) PRINT *, "Non - Cartesian variable block found ", &
              "in file, ignoring and continuing"
          CALL cfd_skip_block()
          CYCLE
        END IF

        ! We now have a valid variable, let's load it up
        ! First error trapping
        CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extent, &
            stagger, mesh_name, mesh_class)

        IF (dims(1) /= nx_global+1 &
            .OR. dims(2) /= ny_global+1 .OR. dims(3) /= nz_global+1) THEN
          IF (rank == 0) PRINT *, "Size of grid represented by one more ", &
              "variables invalid. Continuing"
          CALL cfd_skip_block
          CYCLE
        END IF

        IF (sof /= num) THEN
          IF (rank == 0) PRINT *, "Precision of data does not match ", &
              "precision of code. Continuing."
          CALL cfd_skip_block
        END IF

        ! We're not interested in the other parameters, so if we're here,
        ! load up the data

        CALL cfd_get_3d_cartesian_variable_parallel(data, subtype)

        ! Now have the data, just copy it to correct place

        IF (str_cmp(name(1:3), "Rho")) THEN
          rho(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:6), "Energy")) THEN
          energy(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "Vx")) THEN
          vx(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "Vy")) THEN
          vy(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "Vz")) THEN
          vz(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "Bx")) THEN
          bx(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "By")) THEN
          by(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:2), "Bz")) THEN
          bz(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:6), "p_visc")) THEN
          p_visc(0:nx, 0:ny, 0:nz) = data
        END IF

        IF (str_cmp(name(1:9), "visc_heat")) THEN
          visc_heat(0:nx, 0:ny, 0:nz) = data
        END IF
        
        ! Should be at end of block, but force the point anyway
        CALL cfd_skip_block()
      ELSE
        ! Unknown block, just skip it
        CALL cfd_skip_block()
      END IF
    END DO

    DEALLOCATE(data)

    CALL cfd_close()

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE restart_data

END MODULE setup
