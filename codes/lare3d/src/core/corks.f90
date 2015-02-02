MODULE corks

  USE shared_data

  IMPLICIT NONE
  
  TYPE CORK
    ! id is the 1d (whole) grid index for the cork's starting positon
    INTEGER :: id
    ! the indices of the subgrid cell that contains the cork
    INTEGER :: xi, yi, zi
    ! current position
    REAL(num) :: x, y, z
    ! timestep and time at end of step
    REAL(num) :: dt, t
    ! distance cork has moved over timestep and since deployment
    REAL(num) :: ds, s
    ! is true if cork has just joined subgrid
    LOGICAL :: new
    ! plasma properties are averaged wrt to the cork's position
    REAL(num) :: rho, energy, visc, ohmic, eta
    REAL(num) :: vx, vy, vz
    REAL(num) :: bx, by, bz      
    REAL(num) :: jx, jy, jz    
  END TYPE CORK  

  TYPE CORK_PTR
    TYPE (CORK), POINTER :: crk_p
  END TYPE CORK_PTR


  REAL(num), PARAMETER :: dpy_org_x = 0.0_num
  REAL(num), PARAMETER :: dpy_org_y = 0.25_num
  REAL(num), PARAMETER :: dpy_org_z = 0.0_num  
  REAL(num), PARAMETER :: dpy_rad_min = 0.2_num
  REAL(num), PARAMETER :: dpy_rad_max = 0.4_num
  
  INTEGER, PARAMETER :: cork_fleet_size_max = 50
  INTEGER :: cork_fleet_size, cork_escaped_cnt
  
  TYPE (CORK_PTR), DIMENSION(:), ALLOCATABLE :: cork_fleet, cork_escaped
  
  INTEGER, PARAMETER :: bcast_msg_len = 1
  INTEGER, PARAMETER :: crk_msg_len = 8

  INTEGER, PARAMETER :: crk_df_un = 50
  INTEGER, PARAMETER :: crk_df_field_cnt = 23
  CHARACTER*14, parameter :: crk_df_fn = 'data/corks.dat'
  INTEGER :: crk_df_ios
    

  PRIVATE :: dpy_org_x, dpy_org_y, dpy_org_z
  PRIVATE :: dpy_rad_min, dpy_rad_max
  
  PRIVATE :: cork_fleet_size_max
  PRIVATE :: cork_fleet_size, cork_escaped_cnt
  PRIVATE :: cork_fleet, cork_escaped

  PRIVATE :: bcast_msg_len, crk_msg_len
  PRIVATE :: crk_df_un, crk_df_field_cnt, crk_df_fn, crk_df_ios


CONTAINS
  
  FUNCTION within_grid(x, y, z)
    REAL(num), INTENT(IN) :: x, y, z
    LOGICAL :: within_grid
    within_grid = (x .GE. x_start .AND. x .LE. x_end) .AND. (y .GE. y_start .AND. y .LE. y_end) .AND. (z .GE. z_start .AND. z .LE. z_end)
  END FUNCTION within_grid

  FUNCTION within_grid_ptr(crk_p)
    TYPE (CORK), POINTER :: crk_p
    LOGICAL :: within_grid_ptr

    within_grid_ptr = .FALSE.
    IF (ASSOCIATED(crk_p)) THEN
      within_grid_ptr = within_grid(crk_p%x, crk_p%y, crk_p%z)
    END IF
  END FUNCTION within_grid_ptr

  FUNCTION within_subgrid(x, y, z)
    REAL(num), INTENT(IN) :: x, y, z
    LOGICAL :: within_subgrid

    within_subgrid = (x .GE. xb(0) .AND. x .LE. xb(nx)) .AND. (y .GE. yb(0) .AND. y .LE. yb(ny)) .AND. (z .GE. zb(0) .AND. z .LE. zb(nz))

    ! if point is on the maximum (ni) subgrid boundary (plane, edge or vertex) then within_subgrid is true only if the subgrid boundary
    ! is also the grid boundary
    IF (within_subgrid) THEN
      IF (x .EQ. xb(nx)) THEN
        within_subgrid = x .EQ. x_end
        IF (y .EQ. yb(ny)) THEN
          within_subgrid = within_subgrid .AND. y .EQ. y_end
          IF (z .EQ. zb(nz)) THEN
            within_subgrid = within_subgrid .AND. z .EQ. z_end
          ENDIF
        ELSE IF (z .EQ. zb(nz)) THEN
          within_subgrid = within_subgrid .AND. z .EQ. z_end
        ENDIF
      ELSE IF (y .EQ. yb(ny)) THEN
        within_subgrid = y .EQ. y_end
        IF (z .EQ. zb(nz)) THEN
          within_subgrid = within_subgrid .AND. z .EQ. z_end
        ENDIF
      ELSE IF (z .EQ. zb(nz)) THEN
        within_subgrid = z .EQ. z_end
      END IF 
    END IF

  END FUNCTION within_subgrid

  FUNCTION within_subgrid_ptr(crk_p)
    TYPE (CORK), POINTER :: crk_p
    LOGICAL :: within_subgrid_ptr

    within_subgrid_ptr = .FALSE.
    IF (ASSOCIATED(crk_p)) THEN
      within_subgrid_ptr = within_subgrid(crk_p%x, crk_p%y, crk_p%z)
    END IF
  END FUNCTION within_subgrid_ptr


  ! calculate the global cell index for the specified coordinates
  FUNCTION get_id(x, y, z)
    REAL(num), INTENT(IN) :: x, y, z
    INTEGER :: get_id
    INTEGER :: xi, yi, zi
    REAL(num) :: dx, dy, dz    

    ! calculate cell dimensions
    dx = (x_end - x_start)/nx_global
    dy = (y_end - y_start)/ny_global
    dz = (z_end - z_start)/nz_global

    ! calculate the number of cells between the minimum bound
    ! and the coordinate to the nearest integer
    xi = NINT((x - x_start)/dx)
    yi = NINT((y - y_start)/dy)
    zi = NINT((z - z_start)/dz)    

    ! convert 3d global cell indices into 1d zero-based global index
    get_id = zi*nx_global*ny_global + yi*nx_global + xi   
  END FUNCTION get_id      


  FUNCTION average_over_cube(xfrac,yfrac,zfrac,cube)
    REAL(num), INTENT(IN) :: xfrac, yfrac, zfrac
    REAL(num), DIMENSION(0:1,0:1,0:1), INTENT(IN) :: cube
    REAL(num), DIMENSION(0:1,0:1) :: square
    REAL(num), DIMENSION(0:1) :: line
    REAL(num) :: average_over_cube

    square(0,0) = cube(0,0,0) + (cube(1,0,0)-cube(0,0,0))*xfrac
    square(0,1) = cube(0,0,1) + (cube(1,0,1)-cube(0,0,1))*xfrac
    square(1,1) = cube(0,1,1) + (cube(1,1,1)-cube(0,1,1))*xfrac
    square(1,0) = cube(0,1,0) + (cube(1,1,0)-cube(0,1,0))*xfrac
        
    line(0) = square(0,0) + (square(1,0)-square(0,0))*yfrac
    line(1) = square(0,1) + (square(1,1)-square(0,1))*yfrac 

    average_over_cube = line(0) + (line(1)-line(0))*zfrac
  END FUNCTION average_over_cube

  FUNCTION average_over_vertices(x, y, z, xi, yi, zi, v3d)
    REAL(num), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: xi, yi, zi
    REAL(num), DIMENSION(0:1,0:1,0:1), INTENT(IN) :: v3d
    REAL(num) :: average_over_vertices
    REAL(num) :: xfrac, yfrac, zfrac

    xfrac = ABS(x-xb(xi-1))/ABS(xb(xi)-xb(xi-1))
    yfrac = ABS(y-yb(yi-1))/ABS(yb(yi)-yb(yi-1))
    zfrac = ABS(z-zb(zi-1))/ABS(zb(zi)-zb(zi-1))

    average_over_vertices = average_over_cube(xfrac,yfrac,zfrac,v3d)
  END FUNCTION average_over_vertices

  FUNCTION average_over_face_centres(x, y, z, xi, yi, zi, fi, fc3d)
    REAL(num), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: xi, yi, zi, fi
    REAL(num), DIMENSION(0:1,0:1,0:1), INTENT(IN) :: fc3d
    REAL(num) :: average_over_face_centres
    REAL(num) :: xfrac, yfrac, zfrac
    
    IF (fi .EQ. 0) THEN
      ! x face
      xfrac = ABS(x-xb(xi-1))/ABS(xb(xi)-xb(xi-1))
      yfrac = ABS(y-yc(yi-1))/ABS(yc(yi)-yc(yi-1))
      zfrac = ABS(z-zc(zi-1))/ABS(zc(zi)-zc(zi-1))
    ELSE IF (fi .EQ. 1) THEN
      ! y face
      xfrac = ABS(x-xc(xi-1))/ABS(xc(xi)-xc(xi-1))
      yfrac = ABS(y-yb(yi-1))/ABS(yb(yi)-yb(yi-1))
      zfrac = ABS(z-zc(zi-1))/ABS(zc(zi)-zc(zi-1))
    ELSE
      ! z face
      xfrac = ABS(x-xc(xi-1))/ABS(xc(xi)-xc(xi-1))
      yfrac = ABS(y-yc(yi-1))/ABS(yc(yi)-yc(yi-1))
      zfrac = ABS(z-zb(zi-1))/ABS(zb(zi)-zb(zi-1))
    END IF

    average_over_face_centres = average_over_cube(xfrac,yfrac,zfrac,fc3d)
  END FUNCTION average_over_face_centres

  FUNCTION average_over_cell_centres(x, y, z, xi, yi, zi, cc3d)
    REAL(num), INTENT(IN) :: x, y, z
    INTEGER, INTENT(IN) :: xi, yi, zi
    REAL(num), DIMENSION(0:1,0:1,0:1), INTENT(IN) :: cc3d
    REAL(num) :: average_over_cell_centres
    REAL(num) :: xfrac, yfrac, zfrac

    xfrac = ABS(x-xc(xi-1))/ABS(xc(xi)-xc(xi-1))
    yfrac = ABS(y-yc(yi-1))/ABS(yc(yi)-yc(yi-1))
    zfrac = ABS(z-zc(zi-1))/ABS(zc(zi)-zc(zi-1))
    
    average_over_cell_centres = average_over_cube(xfrac,yfrac,zfrac,cc3d)
  END FUNCTION average_over_cell_centres


  FUNCTION calc_jx(xi, yi, zi)
    INTEGER, INTENT(IN) :: xi, yi, zi
    REAL(num) :: jx1, jx2
    REAL(num) :: calc_jx

    jx1 = (bz(xi,yi+1,zi)-bz(xi,yi,zi))/dyc(yi) - (by(xi,yi,zi+1)-by(xi,yi,zi))/dzc(zi)
    jx2 = (bz(xi+1,yi+1,zi)-bz(xi+1,yi,zi))/dyc(yi) - (by(xi+1,yi,zi+1)-by(xi+1,yi,zi))/dzc(zi)
    
    calc_jx = (jx1 + jx2)*0.5_num
  END FUNCTION calc_jx

  FUNCTION calc_jy(xi, yi, zi)
    INTEGER, INTENT(IN) :: xi, yi, zi
    REAL(num) :: jy1, jy2
    REAL(num) :: calc_jy

    jy1 = (bx(xi,yi,zi+1)-bx(xi,yi,zi))/dzc(zi) - (bz(xi+1,yi,zi)-bz(xi,yi,zi))/dxc(xi)
    jy2 = (bx(xi,yi+1,zi+1)-bx(xi,yi+1,zi))/dzc(zi) - (bz(xi+1,yi+1,zi)-bz(xi,yi+1,zi))/dxc(xi)
    
    calc_jy = (jy1 + jy2)*0.5_num
  END FUNCTION calc_jy

  FUNCTION calc_jz(xi, yi, zi)
    INTEGER, INTENT(IN) :: xi, yi, zi
    REAL(num) :: jz1, jz2
    REAL(num) :: calc_jz

    jz1 = (by(xi+1,yi,zi)-by(xi,yi,zi))/dxc(xi) - (bx(xi,yi+1,zi)-bx(xi,yi,zi))/dyc(yi)
    jz2 = (by(xi+1,yi,zi+1)-by(xi,yi,zi+1))/dxc(xi) - (bx(xi,yi+1,zi+1)-bx(xi,yi,zi+1))/dyc(yi)
    
    calc_jz = (jz1 + jz2)*0.5_num
  END FUNCTION calc_jz


  ! locate the cork within the subgrid
  SUBROUTINE locate_cork(crk_p)
    TYPE (CORK), POINTER, INTENT(IN) :: crk_p
    INTEGER :: xi, yi, zi
    INTEGER :: crk_xi, crk_yi, crk_zi

    IF (ASSOCIATED(crk_p)) THEN

      IF (crk_p%new) THEN
        crk_xi = 1
        crk_yi = 1
        crk_zi = 1
      ELSE
        crk_xi = crk_p%xi
        crk_yi = crk_p%yi
        crk_zi = crk_p%zi
      END IF

      DO xi = 1,nx
        IF (ABS(xc(xi)-crk_p%x) .LT. ABS(xc(crk_xi)-crk_p%x)) THEN
          crk_xi = xi
        END IF
      END DO
      DO yi = 1,ny
        IF (ABS(yc(yi)-crk_p%y) .LT. ABS(yc(crk_yi)-crk_p%y)) THEN
          crk_yi = yi
        END IF
      END DO
      DO zi = 1,nz
        IF (ABS(zc(zi)-crk_p%z) .LT. ABS(zc(crk_zi)-crk_p%z)) THEN
          crk_zi = zi
        END IF
      END DO      

      crk_p%xi = crk_xi
      crk_p%yi = crk_yi
      crk_p%zi = crk_zi  
    
    END IF
  END SUBROUTINE locate_cork


  ! label the cork with the properties of the plasma at its current position
  SUBROUTINE label_cork(crk_p)

    TYPE (CORK), POINTER, INTENT(IN) :: crk_p
    
    INTEGER :: xi, yi, zi
    
    REAL(num), DIMENSION(0:1,0:1,0:1) :: vx3d, vy3d, vz3d
    REAL(num), DIMENSION(0:1,0:1,0:1) :: jx3d, jy3d, jz3d
    REAL(num), DIMENSION(0:1,0:1,0:1) :: bx3d, by3d, bz3d
    REAL(num), DIMENSION(0:1,0:1,0:1) :: rho3d, energy3d
    REAL(num), DIMENSION(0:1,0:1,0:1) :: visc3d, ohmic3d

    IF (ASSOCIATED(crk_p)) THEN
                
      ! find the subgrid cell that contains the cork
      CALL locate_cork(crk_p)
         
      xi = crk_p%xi
      yi = crk_p%yi
      zi = crk_p%zi

      ! average the properties defined at cell vertices for the current cork position      
      vx3d = RESHAPE((/ vx(xi-1,yi-1,zi-1),vx(xi,yi-1,zi-1),vx(xi-1,yi,zi-1),vx(xi,yi,zi-1), &
                        vx(xi-1,yi-1,zi),vx(xi,yi-1,zi),vx(xi-1,yi,zi),vx(xi,yi,zi) /), (/ 2, 2, 2 /))
      vy3d = RESHAPE((/ vy(xi-1,yi-1,zi-1),vy(xi,yi-1,zi-1),vy(xi-1,yi,zi-1),vy(xi,yi,zi-1), &
                        vy(xi-1,yi-1,zi),vy(xi,yi-1,zi),vy(xi-1,yi,zi),vy(xi,yi,zi) /), (/ 2, 2, 2 /))
      vz3d = RESHAPE((/ vz(xi-1,yi-1,zi-1),vz(xi,yi-1,zi-1),vz(xi-1,yi,zi-1),vz(xi,yi,zi-1), &
                        vz(xi-1,yi-1,zi),vz(xi,yi-1,zi),vz(xi-1,yi,zi),vz(xi,yi,zi) /), (/ 2, 2, 2 /))

      crk_p%vx = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,vx3d)
      crk_p%vy = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,vy3d)
      crk_p%vz = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,vz3d)
      
      
      jx3d = RESHAPE((/ calc_jx(xi-1,yi-1,zi-1),calc_jx(xi,yi-1,zi-1),calc_jx(xi-1,yi,zi-1),calc_jx(xi,yi,zi-1), &
                        calc_jx(xi-1,yi-1,zi),calc_jx(xi,yi-1,zi),calc_jx(xi-1,yi,zi),calc_jx(xi,yi,zi) /), (/ 2, 2, 2 /))
      jy3d = RESHAPE((/ calc_jy(xi-1,yi-1,zi-1),calc_jy(xi,yi-1,zi-1),calc_jy(xi-1,yi,zi-1),calc_jy(xi,yi,zi-1), &
                        calc_jy(xi-1,yi-1,zi),calc_jy(xi,yi-1,zi),calc_jy(xi-1,yi,zi),calc_jy(xi,yi,zi) /), (/ 2, 2, 2 /))
      jz3d = RESHAPE((/ calc_jz(xi-1,yi-1,zi-1),calc_jz(xi,yi-1,zi-1),calc_jz(xi-1,yi,zi-1),calc_jz(xi,yi,zi-1), &
                        calc_jz(xi-1,yi-1,zi),calc_jz(xi,yi-1,zi),calc_jz(xi-1,yi,zi),calc_jz(xi,yi,zi) /), (/ 2, 2, 2 /))

      crk_p%jx = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,jx3d)
      crk_p%jy = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,jy3d)
      crk_p%jz = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,jz3d)
                         
      crk_p%eta = average_over_vertices(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,eta)

      ! redfine the cell indices with respect to the cell octant that contains the cork        
      xi = crk_p%xi
      IF (crk_p%x .GT. xc(xi)) THEN          
        xi = crk_p%xi + 1
      END IF
      yi = crk_p%yi
      IF (crk_p%y .GT. yc(yi)) THEN
        yi = crk_p%yi + 1
      END IF
      zi = crk_p%zi
      IF (crk_p%z .GT. zc(zi)) THEN
        zi = crk_p%zi + 1
      END IF

      ! average the properties defined at cell face centres for the current cork position
      bx3d = RESHAPE((/ bx(crk_p%xi-1,yi-1,zi-1),bx(crk_p%xi,yi-1,zi-1),bx(crk_p%xi-1,yi,zi-1),bx(crk_p%xi,yi,zi-1), &
                        bx(crk_p%xi-1,yi-1,zi),bx(crk_p%xi,yi-1,zi),bx(crk_p%xi-1,yi,zi),bx(crk_p%xi,yi,zi) /), (/ 2, 2, 2 /))
      by3d = RESHAPE((/ by(xi-1,crk_p%yi-1,zi-1),by(xi,crk_p%yi-1,zi-1),by(xi-1,crk_p%yi,zi-1),by(xi,crk_p%yi,zi-1), &
                        by(xi-1,crk_p%yi-1,zi),by(xi,crk_p%yi-1,zi),by(xi-1,crk_p%yi,zi),by(xi,crk_p%yi,zi) /), (/ 2, 2, 2 /))
      bz3d = RESHAPE((/ bz(xi-1,yi-1,crk_p%zi-1),bz(xi,yi-1,crk_p%zi-1),bz(xi-1,yi,crk_p%zi-1),bz(xi,yi,crk_p%zi-1), &
                        bz(xi-1,yi-1,crk_p%zi),bz(xi,yi-1,crk_p%zi),bz(xi-1,yi,crk_p%zi),bz(xi,yi,crk_p%zi) /), (/ 2, 2, 2 /))

      crk_p%bx = average_over_face_centres(crk_p%x,crk_p%y,crk_p%z,crk_p%xi,yi,zi,0,bx3d)
      crk_p%by = average_over_face_centres(crk_p%x,crk_p%y,crk_p%z,xi,crk_p%yi,zi,1,by3d)
      crk_p%bz = average_over_face_centres(crk_p%x,crk_p%y,crk_p%z,xi,yi,crk_p%zi,2,bz3d)
        

      ! average the properties defined at cell centres for the current cork position
      rho3d = RESHAPE((/ rho(xi-1,yi-1,zi-1),rho(xi,yi-1,zi-1),rho(xi-1,yi,zi-1),rho(xi,yi,zi-1), &
                         rho(xi-1,yi-1,zi),rho(xi,yi-1,zi),rho(xi-1,yi,zi),rho(xi,yi,zi) /), (/ 2, 2, 2 /))
      energy3d = RESHAPE((/ energy(xi-1,yi-1,zi-1),energy(xi,yi-1,zi-1),energy(xi-1,yi,zi-1),energy(xi,yi,zi-1), &
                            energy(xi-1,yi-1,zi),energy(xi,yi-1,zi),energy(xi-1,yi,zi),energy(xi,yi,zi) /), (/ 2, 2, 2 /))
      visc3d = RESHAPE((/ visc_heat(xi-1,yi-1,zi-1),visc_heat(xi,yi-1,zi-1),visc_heat(xi-1,yi,zi-1),visc_heat(xi,yi,zi-1), &
                          visc_heat(xi-1,yi-1,zi),visc_heat(xi,yi-1,zi),visc_heat(xi-1,yi,zi),visc_heat(xi,yi,zi) /), (/ 2, 2, 2 /))
      ohmic3d = RESHAPE((/ ohmic_heat(xi-1,yi-1,zi-1),ohmic_heat(xi,yi-1,zi-1),ohmic_heat(xi-1,yi,zi-1),ohmic_heat(xi,yi,zi-1), &
                           ohmic_heat(xi-1,yi-1,zi),ohmic_heat(xi,yi-1,zi),ohmic_heat(xi-1,yi,zi),ohmic_heat(xi,yi,zi) /), (/ 2, 2, 2 /))

      crk_p%rho = average_over_cell_centres(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,rho3d) 
      crk_p%energy = average_over_cell_centres(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,energy3d) 
      crk_p%visc = average_over_cell_centres(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,visc3d) 
      crk_p%ohmic = average_over_cell_centres(crk_p%x,crk_p%y,crk_p%z,xi,yi,zi,ohmic3d) 
      
    END IF ! ASSOCIATED(crk_p)

  END SUBROUTINE label_cork


  ! cork is deployed (starting position must be within subgrid) at the position given by the input parameters
  ! plasma properties are recorded
  SUBROUTINE deploy_cork(x, y, z)
    REAL(num), INTENT(IN) :: x, y, z
    
    TYPE (CORK), POINTER :: crk_p
      
    IF (cork_fleet_size .LT. cork_fleet_size_max) THEN  
      IF (within_subgrid(x,y,z)) THEN 
                  
        ALLOCATE(crk_p)      

        IF (ASSOCIATED(crk_p)) THEN
          
          crk_p%id = get_id(x,y,z)
          crk_p%x = x
          crk_p%y = y
          crk_p%z = z
          crk_p%dt = dt
          crk_p%t = time
          crk_p%ds = 0.0_num          
          crk_p%s = 0.0_num

          crk_p%new = .TRUE.

          CALL label_cork(crk_p)
                    
          ! add the cork to the fleet
          cork_fleet_size = cork_fleet_size + 1    
          cork_fleet(cork_fleet_size)%crk_p => crk_p

          PRINT *, "Cork ", crk_p%id, " deployed by process ", rank, "."
          
        END IF ! ASSOCIATED(crk_p)

      END IF ! within_subgrid(x,y,z)      
    END IF ! cork_fleet_size .LT. cork_fleet_size_max

  END SUBROUTINE deploy_cork

  
  ! deploy corks at random locations within a spherical shell centred at dpy_org_xyz and of thickness dpy_rad_max-dpy_rad_min 
  SUBROUTINE deploy_cork_fleet
    INTEGER :: crk_i, dpy_i, crd_i, id, rank_zero

    REAL(num), DIMENSION(1:3) :: dpy_rnd    
    REAL(num) :: rad, theta, phi
    REAL(num), DIMENSION(:), ALLOCATABLE :: dpy_pos
    INTEGER, DIMENSION(:), ALLOCATABLE :: dpy_id
           
    ALLOCATE(cork_fleet(1:cork_fleet_size_max))    
    ALLOCATE(cork_escaped(1:cork_fleet_size_max))        
    DO crk_i = 1,cork_fleet_size_max
      NULLIFY(cork_fleet(crk_i)%crk_p)
      NULLIFY(cork_escaped(crk_i)%crk_p)
    END DO
    cork_fleet_size = 0
    cork_escaped_cnt = 0      

    ALLOCATE(dpy_pos(1:cork_fleet_size_max*3))
    rank_zero = 0

    IF (rank .EQ. rank_zero) THEN
      ! create corks.dat file    
      OPEN(UNIT=crk_df_un, STATUS='NEW', FILE=crk_df_fn, FORM='UNFORMATTED', ACCESS='STREAM', iostat=crk_df_ios)
      WRITE(crk_df_un) num, crk_df_field_cnt 
      CLOSE(UNIT=crk_df_un)
                
      ALLOCATE(dpy_id(1:cork_fleet_size_max))
      CALL RANDOM_SEED      
      DO dpy_i = 1,cork_fleet_size_max
        crd_i = (dpy_i-1)*3
 
        dpy_id(dpy_i) = -1
        DO WHILE (dpy_id(dpy_i) .EQ. -1)
          CALL RANDOM_NUMBER(dpy_rnd)

          rad = dpy_rad_min + dpy_rnd(1)*(dpy_rad_max-dpy_rad_min)
          theta = dpy_rnd(2)*pi
          phi = dpy_rnd(3)*2.0_num*pi

          dpy_pos(crd_i+1) = rad*SIN(theta)*COS(phi) + dpy_org_x
          dpy_pos(crd_i+2) = rad*SIN(theta)*SIN(phi) + dpy_org_y
          dpy_pos(crd_i+3) = rad*COS(theta) + dpy_org_z

          id = get_id(dpy_pos(crd_i+1), dpy_pos(crd_i+2), dpy_pos(crd_i+3))
          IF (ANY(dpy_id(1:dpy_i) .NE. id)) THEN
            ! each cork must be deployed to a unique grid cell
            dpy_id(dpy_i) = id
          END IF
        END DO         
      END DO      
      DEALLOCATE(dpy_id)

      CALL MPI_BCAST(dpy_pos, cork_fleet_size_max*3, MPI_DOUBLE_PRECISION, rank_zero, comm, errcode)
      IF (errcode .NE. MPI_SUCCESS) THEN
        PRINT *, "Error ", errcode, " broadcasting deployment positions from rank ", rank_zero, "."
      END IF  
    ELSE
      CALL MPI_BCAST(dpy_pos, cork_fleet_size_max*3, MPI_DOUBLE_PRECISION, rank_zero, comm, errcode)
      IF (errcode .NE. MPI_SUCCESS) THEN
        PRINT *, "Error ", errcode, " receiving deployment positions from rank ", rank_zero, " to rank ", rank, "."
      END IF  
    END IF
    
    DO dpy_i = 1,cork_fleet_size_max      
      crd_i = (dpy_i-1)*3
      CALL deploy_cork(dpy_pos(crd_i+1), dpy_pos(crd_i+2), dpy_pos(crd_i+3))      
    END DO  

    DEALLOCATE(dpy_pos)

    CALL MPI_BARRIER(comm, errcode)

    CALL output_cork_fleet
    
  END SUBROUTINE deploy_cork_fleet

  
  ! dellocate the memory used to store the cork fleet
  SUBROUTINE retire_cork_fleet
    INTEGER :: crk_i
    
    DO crk_i = 1,cork_fleet_size
      IF (ASSOCIATED(cork_fleet(crk_i)%crk_p)) THEN
        DEALLOCATE(cork_fleet(crk_i)%crk_p)
        NULLIFY(cork_fleet(crk_i)%crk_p)
      END IF
      IF (ASSOCIATED(cork_escaped(crk_i)%crk_p)) THEN
        DEALLOCATE(cork_escaped(crk_i)%crk_p)
        NULLIFY(cork_escaped(crk_i)%crk_p)
      END IF      
    END DO    

    DEALLOCATE(cork_fleet)    
    DEALLOCATE(cork_escaped)    
    
    cork_fleet_size = 0
    cork_escaped_cnt = 0

  END SUBROUTINE retire_cork_fleet
   
 
  ! add a cork to the fleet
  SUBROUTINE add_cork(crk_msg)
    REAL(num), DIMENSION(1:crk_msg_len), INTENT(IN) :: crk_msg

    INTEGER :: crk_i
    TYPE (CORK), POINTER :: crk_p
      
    IF (cork_fleet_size .LT. cork_fleet_size_max) THEN
      crk_i = cork_fleet_size + 1      
      ALLOCATE(cork_fleet(crk_i)%crk_p)
      crk_p => cork_fleet(crk_i)%crk_p
      IF (ASSOCIATED(crk_p)) THEN
        crk_p%id = INT(crk_msg(1))
        crk_p%x = crk_msg(2)
        crk_p%y = crk_msg(3)
        crk_p%z = crk_msg(4)

        crk_p%dt = crk_msg(5)
        crk_p%t = crk_msg(6)
        crk_p%ds = crk_msg(7)
        crk_p%s = crk_msg(8)
        
        crk_p%new = .TRUE.

        CALL label_cork(crk_p)
              
        cork_fleet_size = cork_fleet_size + 1
      END IF
    END IF    

  END SUBROUTINE add_cork

    
  ! every processor broadcasts the corks that have left its subgrid (but not the grid)
  ! every processor receives the corks that have entered its subgrid
  SUBROUTINE reform_cork_fleet

    INTEGER :: crk_i
    TYPE (CORK), POINTER :: crk_p
    
    INTEGER :: broadcaster    
    INTEGER, DIMENSION(1:bcast_msg_len) :: bcast_msg
    REAL(num), DIMENSION(1:crk_msg_len) :: crk_msg
        

    DO broadcaster = 0,nproc-1
      IF (rank .EQ. broadcaster) THEN

        bcast_msg = (/ cork_escaped_cnt /)
        CALL MPI_BCAST(bcast_msg, bcast_msg_len, MPI_INTEGER, broadcaster, comm, errcode)
        IF (errcode .NE. MPI_SUCCESS) THEN
          PRINT *, "Error ", errcode, " broadcasting escaped cork cnt from ", rank, "."
        ELSE
          DO crk_i = 1,cork_escaped_cnt 
            crk_p => cork_escaped(crk_i)%crk_p
            IF (ASSOCIATED(crk_p)) THEN
              !PRINT *, "Process ", rank, ": broadcasting escaped cork ", crk_p%id, "."
              crk_msg = (/ DBLE(crk_p%id), crk_p%x, crk_p%y, crk_p%z, crk_p%dt, crk_p%t, crk_p%ds, crk_p%s /)
              CALL MPI_BCAST(crk_msg, crk_msg_len, MPI_DOUBLE_PRECISION, broadcaster, comm, errcode)
              IF (errcode .NE. MPI_SUCCESS) THEN
                PRINT *, "Error ", errcode, " broadcasting escaped cork position from ", rank, "."                
              END IF
            END IF
          END DO
        END IF

      ELSE

        CALL MPI_BCAST(bcast_msg, bcast_msg_len, MPI_INTEGER, broadcaster, comm, errcode)
        IF (errcode .NE. MPI_SUCCESS) THEN
          PRINT *, "Error ", errcode, " receiving escaped cork cnt from ", broadcaster, " to ", rank, "."
          bcast_msg(1) = 0
        ELSE
          !PRINT *, "Process ", rank, ": receiving ", bcast_msg(1), " escaped cork(s)."
          DO crk_i = 1,bcast_msg(1)            
            CALL MPI_BCAST(crk_msg, crk_msg_len, MPI_DOUBLE_PRECISION, broadcaster, comm, errcode)
            IF (errcode .NE. MPI_SUCCESS) THEN
              PRINT *, "Error ", errcode, " receiving new cork position from ", broadcaster, " to ", rank, "."
            ELSE
              !PRINT *, "Process ", rank, ": receiving escaped cork ", INT(crk_msg(1)), "."
              IF (within_subgrid(crk_msg(2),crk_msg(3),crk_msg(4))) THEN
                ! cork has moved into subgrid so add it to fleet                
                !PRINT *, "Process ", rank, ": added escaped cork ", INT(crk_msg(1)), "."
                CALL add_cork(crk_msg)
              END IF              
            END IF
          END DO
        END IF

      END IF ! rank .EQ. broadcaster

      CALL MPI_BARRIER(comm, errcode)

    END DO ! DO broadcaster = 0,nproc-1

  END SUBROUTINE reform_cork_fleet


  ! ensure the next free space in the escaped array points to the escaped cork
  SUBROUTINE escape_cork(crk_i)
    INTEGER, INTENT(IN) :: crk_i

    TYPE (CORK), POINTER :: crk_p

    IF (crk_i .GE. 1 .AND. crk_i .LE. cork_fleet_size) THEN
      IF (ALLOCATED(cork_escaped)) THEN
        IF (cork_escaped_cnt .LT. SIZE(cork_escaped,1)) THEN
          crk_p => cork_fleet(crk_i)%crk_p
          IF (ASSOCIATED(crk_p)) THEN
            cork_escaped_cnt = cork_escaped_cnt + 1
            NULLIFY(cork_escaped(cork_escaped_cnt)%crk_p)
            cork_escaped(cork_escaped_cnt)%crk_p => crk_p
          END IF
        END IF
      END IF
    END IF

  END SUBROUTINE escape_cork


  ! update the cork's position using the velocities defined at the half time step
  ! if cork moves out of subgrid broadcast new position and remove cork from fleet
  ! add corks to fleet that have moved into subgrid
  ! record plasma properties at new cork position(s) 
  SUBROUTINE update_cork_fleet()

    INTEGER :: crk_i, crk_i2
    INTEGER :: xi, yi, zi
           
    TYPE (CORK), POINTER :: crk_p
            
    REAL(num), DIMENSION(0:1,0:1,0:1) :: vx3d, vy3d, vz3d
    REAL(num) :: x, y, z
    
        
    ! update cork positions and remove any corks that have left the subgrid    
    DO crk_i = 1,cork_fleet_size

      crk_p => cork_fleet(crk_i)%crk_p

      IF (ASSOCIATED(crk_p)) THEN
                
        !PRINT *, "Process ", rank, " updating cork ", crk_p%id, "."

        ! the cell that contains the cork
        xi = crk_p%xi
        yi = crk_p%yi
        zi = crk_p%zi

        vx3d = RESHAPE((/ vx1(xi-1,yi-1,zi-1),vx1(xi,yi-1,zi-1),vx1(xi-1,yi,zi-1),vx1(xi,yi,zi-1), &
                          vx1(xi-1,yi-1,zi),vx1(xi,yi-1,zi),vx1(xi-1,yi,zi),vx1(xi,yi,zi) /), (/ 2, 2, 2 /))
        vy3d = RESHAPE((/ vy1(xi-1,yi-1,zi-1),vy1(xi,yi-1,zi-1),vy1(xi-1,yi,zi-1),vy1(xi,yi,zi-1), &
                          vy1(xi-1,yi-1,zi),vy1(xi,yi-1,zi),vy1(xi-1,yi,zi),vy1(xi,yi,zi) /), (/ 2, 2, 2 /))
        vz3d = RESHAPE((/ vz1(xi-1,yi-1,zi-1),vz1(xi,yi-1,zi-1),vz1(xi-1,yi,zi-1),vz1(xi,yi,zi-1), &
                          vz1(xi-1,yi-1,zi),vz1(xi,yi-1,zi),vz1(xi-1,yi,zi),vz1(xi,yi,zi) /), (/ 2, 2, 2 /))

        ! the current cell position
        x = crk_p%x
        y = crk_p%y
        z = crk_p%z

        crk_p%x = x + average_over_vertices(x,y,z,xi,yi,zi,vx3d)*dt
        crk_p%y = y + average_over_vertices(x,y,z,xi,yi,zi,vy3d)*dt
        crk_p%z = z + average_over_vertices(x,y,z,xi,yi,zi,vz3d)*dt

        crk_p%dt = dt
        crk_p%t = time
        crk_p%ds = SQRT((crk_p%x-x)**2 + (crk_p%y-y)**2 + (crk_p%z-z)**2)
        crk_p%s = crk_p%s + crk_p%ds

        crk_p%new = .FALSE.
         
        IF (within_subgrid_ptr(crk_p)) THEN
          !PRINT *, "Process ", rank, ": cork ", crk_p%id, " within subgrid."
          CALL label_cork(crk_p)
        ELSE
          IF (within_grid_ptr(crk_p)) THEN
            !PRINT *, "Process ", rank, ": cork ", crk_p%id, " has escaped."
            CALL escape_cork(crk_i)
          ELSE
            PRINT *, "Process ", rank, ": cork ", crk_p%id, " is lost."          
          END IF
          NULLIFY(cork_fleet(crk_i)%crk_p)
        END IF
        
      END IF
    END DO
    
    ! remove the spaces in the fleet that have been left by escaped and lost corks
    cork_fleet_size = 0
    DO crk_i = 1,cork_fleet_size_max
      IF (ASSOCIATED(cork_fleet(crk_i)%crk_p)) THEN
        cork_fleet_size = cork_fleet_size + 1
      ELSE
        DO crk_i2 = crk_i,cork_fleet_size_max-1
          cork_fleet(crk_i2)%crk_p => cork_fleet(crk_i2+1)%crk_p
        END DO
        NULLIFY(cork_fleet(cork_fleet_size_max)%crk_p)        
      END IF
    END DO

    CALL MPI_BARRIER(comm, errcode)
    ! adjust fleet in case any corks have left or entered the subgrid
    CALL reform_cork_fleet    
    CALL MPI_BARRIER(comm, errcode)

    ! deallocate memory used for escaped corks
    DO crk_i = 1,cork_escaped_cnt
      IF (ASSOCIATED(cork_escaped(crk_i)%crk_p)) THEN
        DEALLOCATE(cork_escaped(crk_i)%crk_p)
        NULLIFY(cork_escaped(crk_i)%crk_p)
      END IF
    END DO
    cork_escaped_cnt = 0       

    CALL output_cork_fleet
       
  END SUBROUTINE update_cork_fleet


  ! write the current cork data (e.g, position and local plasma properties) to a data file
  SUBROUTINE output_cork(crk_p)

    TYPE (CORK), POINTER, INTENT(IN) :: crk_p
    INTEGER :: crk_df_pos
    REAL(num) :: dble_new

    IF (ASSOCIATED(crk_p)) THEN
                        
      dble_new = 0.0_num
      IF (crk_p%new) THEN
        dble_new = 1.0_num
      END IF

      OPEN(UNIT=crk_df_un, STATUS='OLD', FILE=crk_df_fn, FORM='UNFORMATTED', ACCESS='STREAM', ACTION='WRITE', iostat=crk_df_ios)
      INQUIRE(UNIT=crk_df_un, SIZE=crk_df_pos)
      WRITE(crk_df_un, POS=crk_df_pos+1) DBLE(crk_p%id), dble_new
      WRITE(crk_df_un) crk_p%x, crk_p%y, crk_p%z
      WRITE(crk_df_un) crk_p%dt, crk_p%t, crk_p%ds, crk_p%s 
      WRITE(crk_df_un) crk_p%rho, crk_p%energy, crk_p%visc, crk_p%ohmic, crk_p%eta
      WRITE(crk_df_un) crk_p%vx, crk_p%vy, crk_p%vz, crk_p%bx, crk_p%by, crk_p%bz, crk_p%jx, crk_p%jy, crk_p%jz
      CLOSE(UNIT=crk_df_un)
     
    END IF 
        
  END SUBROUTINE output_cork


  SUBROUTINE output_cork_fleet
    INTEGER :: writer, crk_i
    TYPE (CORK), POINTER :: crk_p

    DO writer = 0,nproc-1
      IF (rank .EQ. writer) THEN
        DO crk_i = 1,cork_fleet_size
          crk_p => cork_fleet(crk_i)%crk_p
          IF (ASSOCIATED(crk_p)) THEN
            CALL output_cork(crk_p)
          END IF
        END DO
      END IF

      CALL MPI_BARRIER(comm, errcode)
    END DO
    
  END SUBROUTINE output_cork_fleet
  


END MODULE corks

  

