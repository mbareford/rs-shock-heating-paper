!****************************************************************
! All global variables defined here (cf F77 COMMON block).
! All the names in here are public provided the MODULE is USE'd
!****************************************************************
MODULE constants

  IMPLICIT NONE

#ifdef Q_SINGLE
  INTEGER, PARAMETER :: num = KIND(1.0) 
#else
  INTEGER, PARAMETER :: num = KIND(1.D0) 
#endif
  INTEGER, PARAMETER :: dbl = KIND(1.D0)

  REAL(num), PARAMETER :: pi = 3.14159265358979323_num

  ! These are the real SI physical constants
  ! Permiability of free space
  REAL(num), PARAMETER :: mu0_si =  4.0e-7_num * pi
  
  ! Gas Constant
  REAL(num), PARAMETER :: gc_si = 8.3144621_num

  ! Avogadro Constant
  REAL(num), PARAMETER :: av_si = 6.02214129e23_num
  
  ! Electron Charge
  REAL(num), PARAMETER :: e_si = 1.60217657e-19_num
  REAL(num), PARAMETER :: e_cgs = 4.8032068e-10_num
  
  ! Permittivity of Free Space
  REAL(num), PARAMETER :: e0_si = 8.8541878e-12_num

  ! Boltzmann's Constant
  REAL(num), PARAMETER :: kb_si = 1.3806504e-23_num
  REAL(num), PARAMETER :: kb_cgs = 1.6e-12_num
  
  ! Mass of proton
  REAL(num), PARAMETER :: mp_si = 1.67262158e-27_num
  
  ! Mass of electron
  REAL(num), PARAMETER :: me_si = 9.10938188e-31_num
  
  ! Planck's constant
  REAL(num), PARAMETER :: hp_si = 6.626068e-34_num
  
  ! Ionisation potential of hydrogen in J
  REAL(num), PARAMETER :: ionise_pot_si = 2.17870364e-18_num

  ! Solar radius
  REAL(num), PARAMETER :: solar_rad_si = 6.955e8_num

  ! Solar mass
  REAL(num), PARAMETER :: solar_mass_si = 1.9891e30_num

  ! Gravitational Constant
  REAL(num), PARAMETER :: gvc_si = 6.67398e-11_num
    

  REAL(num), PARAMETER :: none_zero = TINY(1.0_num) 
  REAL(num), PARAMETER :: largest_number = HUGE(1.0_num)  
  REAL(num), PARAMETER :: third = 1.0_num / 3.0_num, sixth = 1.0_num / 6.0_num
  INTEGER, PARAMETER :: BC_PERIODIC = 1, BC_OTHER = 2
  INTEGER, PARAMETER :: BC_OPEN = 3

  INTEGER, PARAMETER :: version = 2, revision = 10

  ! IC codes
  ! This is a bitmask, remember that
  INTEGER, PARAMETER :: IC_NEW = 1, IC_RESTART = 2

  ! Equation of state codes
  INTEGER, PARAMETER :: EOS_IDEAL = 1, EOS_ION = 2, EOS_PI = 3


END MODULE constants


MODULE shared_data

  USE constants
  IMPLICIT NONE
  INCLUDE 'mpif.h'

#ifdef Q_SINGLE
  INTEGER :: mpireal = MPI_REAL
#else
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
#endif  

  INTEGER :: nx_global, ny_global, nz_global
  ! NB: as there are now 2 ghost celss so indexing will fail if (nx, ny, nz)<2
  INTEGER :: nx, ny, nz
  INTEGER  :: nsteps
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: rho, energy
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bx, vx, vx1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: by, vy, vy1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: bz, vz, vz1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: jx_r, jy_r, jz_r, rho_r
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: fx_r, fy_r, fz_r
  
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: delta_ke, p_visc, visc_heat, ohmic_heat
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: visc_heat_xy, visc_heat_xz, visc_heat_yz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: visc_heat_xx, visc_heat_yy, visc_heat_zz
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eta, eta_num, j_crit, visc_r
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: cv, cv1
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: dis_dzbx_w, dis_dzby_w

  REAL(num), DIMENSION(:), ALLOCATABLE :: xc, xb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dxb, dxc
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc, yb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dyb, dyc, grav
  REAL(num), DIMENSION(:), ALLOCATABLE :: zc, zb
  REAL(num), DIMENSION(:), ALLOCATABLE :: dzb, dzc

  INTEGER, PARAMETER :: data_dir_max_length = 64
  CHARACTER(LEN = data_dir_max_length) :: data_dir

  REAL(num) :: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12
  REAL(num) :: dt, dt2, dtr, dth, t_end, time, dt_multiplier
  REAL(num) :: length_x, length_y, length_z, visc1, visc2, visc3
  REAL(num) :: x_start, x_end, y_start, y_end, z_start, z_end
  REAL(num) :: gamma, eta_crit, dt_snapshots, eta_background, j_crit_const
  REAL(num) :: heating_visc = 0.0_num, heating_ohmic = 0.0_num, heating_dp = 0.0_num
  REAL(num) :: htvisc_xy = 0.0_num, htvisc_xz = 0.0_num, htvisc_yz = 0.0_num
  REAL(num) :: htvisc_xx = 0.0_num, htvisc_yy = 0.0_num, htvisc_zz = 0.0_num
  REAL(num) :: con_supp = 0.0_num, x_loss_con = 0.0_num, y_loss_con = 0.0_num, z_loss_con = 0.0_num, loss_rad = 0.0_num
  REAL(num) :: delta_helicity = 0.0_num, eta_crit_cnt = 0.0_num
  REAL(num) :: tf1=0.0_num,tf2=0.0_num,tf3=0.0_num
  REAL(num) :: tdbe1=0.0_num, tdbe2=0.0_num, tdbe3=0.0_num, tdbe4=0.0_num, tdke1=0.0_num
  REAL(num) :: rke_neg=0.0_num, rke_pos=0.0_num, rbe_neg=0.0_num, rbe_pos=0.0_num
  REAL(num) :: dxby1_w=0.0_num, dxbz1_w=0.0_num, dybx1_w=0.0_num, dybz1_w=0.0_num, dzbx1_w=0.0_num, dzby1_w=0.0_num
  REAL(num) :: dxby2_w=0.0_num, dxbz2_w=0.0_num, dybx2_w=0.0_num, dybz2_w=0.0_num, dzbx2_w=0.0_num, dzby2_w=0.0_num
      
  INTEGER :: xbc_max, ybc_max, xbc_min, ybc_min, zbc_min, zbc_max
  INTEGER :: ix, iy, iz, ixp, iyp, izp, ixm, iym, izm, xpass, ypass, zpass
  INTEGER :: restart_snapshot
  INTEGER :: peak_substeps = 0
  LOGICAL :: x_stretch, y_stretch, z_stretch, rke
  LOGICAL :: resistive_mhd, any_open
  LOGICAL :: restart

  ! normalising constants                 
  REAL(num) :: B0, L0, NmDn0, m0, Dn0, P0, Tp0, vA, En0, tm0, K0, eta0
  ! mass fraction - mass of ions in units of proton mass
  REAL(num) :: mf

  ! initial state  
  REAL(num) :: iniTp, iniEn
  REAL(num) :: lp_rad, lp_alp, lp_omg

  ! diagnostic lagrangian tracers
  LOGICAL :: diag_corks
  
  ! kink perturbation
  LOGICAL :: perturbation
  REAL(num) :: ptb_sign, ptb_k, ptb_amp, ptb_lim 

  ! Conduction
  LOGICAL :: conduction, con_flux_limited
  REAL(num) :: con_limit, con_b_min

  ! Radiation
  LOGICAL :: radiation, rad_implicit, rad_avTp
  INTEGER :: rad_alg
  REAL(num) :: rad_min_tp, rad_min_loss, rad_dTp_frac

  ! Successive over-relaxation
  INTEGER :: sor_it_lim, sor_err_wgt
  REAL(num) :: sor_w, sor_err

  LOGICAL :: conrad_half_dt

  ! Equation of state
  INTEGER :: eos_number = EOS_IDEAL

  ! Damping boundary variables
  LOGICAL :: damping

  ! Partially ionised plasma
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: eta_perp, xi_n, eta_perp0
  REAL(num), DIMENSION(:, :, :), ALLOCATABLE :: parallel_current, perp_current
  LOGICAL :: cowling_resistivity, neutral_gas
  REAL(num) :: fbar, tbar, tr, ionise_pot, rbar
  REAL(num) :: etabar

  ! MPI data
  INTEGER :: rank, proc_x_min, proc_x_max, proc_y_min, proc_y_max, proc_z_min, proc_z_max, coordinates(3)
  INTEGER :: errcode, comm, tag, nproc, nprocx, nprocy, nprocz
  INTEGER :: status(MPI_STATUS_SIZE)

  ! file handling
  INTEGER :: subtype, obstype
  INTEGER(KIND = MPI_OFFSET_KIND) :: initialdisp
  INTEGER :: initial
  INTEGER :: n_zeros = 4
  INTEGER :: output_file = 0

  ! Number of variables to dump
  LOGICAL, DIMENSION(40) :: dump_mask

END MODULE shared_data



! The pre-processor removes the following line so it compiles without error
! unless the pre-processor hasn't been run over it

#ifdef PRECHECK
This line deliberately breaks the compile IF the preprocessor has not worked.
#endif
