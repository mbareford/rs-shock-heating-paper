MODULE control

  USE shared_data
  USE normalise

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: user_normalisation, control_variables, set_output_dumps

CONTAINS

  SUBROUTINE user_normalisation

    ! Gamma is the ratio of specific heat capacities
    gamma = (5.0_num / 3.0_num)

    ! Average mass of an ion in proton masses
    ! The code assumes a single ion species with this mass
    mf = 1.0_num
    m0 = mf*mp_si

    ! The equations describing the normalisation in LARE
    ! have three free parameters which must be specified by
    ! the end user. 

    ! magnetic field
    B0 = 5.0e-3_num
    ! length
    L0 = 1.0e6_num
    ! particle number density
    NmDn0 = 1.0e15_num

    ! Other parameters are also specified.

    ! density
    Dn0 = m0*NmDn0
    ! velocity
    vA = B0/SQRT(mu0_si*Dn0)
    ! time
    tm0 = L0/vA
    ! pressure
    P0 = (B0**2)/mu0_si
    ! temperature
    Tp0 = (P0/Dn0)*(m0/kb_si)
    ! thermal energy
    En0 = P0/Dn0

    ! conduction kappa    
    K0 = 1.e-11_num*(Tp0**3.5_num)*tm0/(Dn0*(L0**2)*En0)

    ! resistivity
    eta0 = (mu0_si*L0*vA)*(Dn0**2)*SQRT(Tp0)/(B0**2)
        
  END SUBROUTINE user_normalisation



  SUBROUTINE control_variables

    !REAL(num) :: r_lp
    
    ! The code to choose the initial conditions. The valid choices are
    ! IC_NEW - Use set_initial_conditions in "initial_conditions.f90"
    !         to setup new initial conditions
    ! IC_RESTART - Load the output file with index restart_snapshot and
    ! use it as the initial conditions
    initial = IC_NEW
    restart_snapshot = 0


    ! Set the number of gridpoints in x and y directions
    nx_global = 256
    ny_global = 256
    nz_global = 512

    ! Set the maximum number of iterations of the core solver before the code
    ! terminates. If nsteps < 0 then the code will run until t = t_end
    nsteps = -1

    ! The maximum runtime of the code
    t_end = 400.0_num

    ! CFL multiplier
    dt_multiplier = 1.0_num/SQRT(3.0_num)


    ! Set these constants to manually
    ! override the domain decomposition.
    ! If either constant is set to zero
    ! then the code will try to automatically
    ! decompose in this direction
    nprocx = 0
    nprocy = 0
    nprocz = 0

    ! The length of the domain in the x direction
    x_start = -2.0_num
    x_end = 2.0_num
    ! Should the x grid be stretched or uniform
    x_stretch = .FALSE.

    y_start = -2.0_num
    y_end = 2.0_num
    y_stretch = .FALSE.

    z_start = -10.0_num
    z_end = 10.0_num
    z_stretch = .FALSE.

    ! Set the boundary conditions on the four edges of the simulation domain
    ! Valid constants are
    ! BC_PERIODIC - Periodic boundary conditions
    ! BC_OPEN - Reimann characteristic boundary conditions
    ! BC_OTHER - User boundary conditions specified in "boundary.f90"
    xbc_min = BC_OTHER
    xbc_max = BC_OTHER
    ybc_max = BC_OTHER
    ybc_min = BC_OTHER
    zbc_min = BC_OTHER
    zbc_max = BC_OTHER


    ! Remap kinetic energy correction. LARE does not
    ! perfectly conserve kinetic energy during the remap step
    ! This missing energy can be added back into the simulation
    ! as a uniform heating. Turning rke to true turns on this
    ! addition
    rke = .TRUE.

    ! If cowling_resistivity is true then the code calculates and
    ! applies the Cowling Resistivity to the MHD equations     
    ! only possible if not EOS_IDEAL
    cowling_resistivity = .FALSE.

    ! set to true to turn on routine for damped boundaries
    ! these routines are in boundary.f90 and you should check these
    ! actually do what you want.
    damping = .FALSE.

    ! Set the equation of state. Valid choices are
    ! EOS_IDEAL - Simple ideal gas for perfectly ionised plasma
    ! EOS_PI - Simple ideal gas for partially ionised plasma
    ! EOS_ION - EOS_PI plus the ionisation potential 
    ! N.B. read the manual for notes on these choices
    eos_number = EOS_IDEAL
    ! EOS_IDEAL also requires that you specific whether
    ! the gas is ionised or not. Some stratified atmospheres
    ! only work for neutral hydrogen for example
    ! For fully ionised gas set .FALSE.
    ! For neutral hydrogen set .TRUE.
    ! This flag is ignored for all other EOS choices.
    neutral_gas = .FALSE.


    ! Shock viscosities as detailed in manual - they are dimensionless
    visc1 = 0.1_num
    visc2 = 0.5_num
    ! Real viscosity expressed as the inverse Reynolds number
    visc3 = 0.0_num


    ! Turn on or off the resistive parts of the MHD equations
    resistive_mhd = .TRUE.

    ! The background resistivity expressed as the inverse Lundquist number
    eta_background = 0.0_num

    ! The anomalous resistivity is expressed as the inverse Lundquist number
    eta_crit = 0.001_num


    ! This constant is used to determine the critical threshold for anomalous resistivity.
    ! It is calculated by equating the electron drift speed with the ion sound speed
    ! and by using the Larmor proton radius to adjust for under-resolved current.
    
    ! the critical current constant is calculated used a local Larmor proton radius
    j_crit_const = (SQRT(gamma)/2.0_num)*((mu0_si*kB_si)/me_si)*((Dn0*Tp0)/(B0**2))

    ! proton larmor radius calculated from a temperature of 1MK    
    !r_lp = SQRT(kb_si*mp_si*1.0e6_num)/(e_si*B0)
    ! the critical current constant is calculated used a fixed Larmor proton radius
    !j_crit_const = 0.5_num*((e_si*mu0_si)/me_si)*(Dn0/B0)*SQRT((gamma*kb_si*Tp0)/mp_si)*r_lp
    
    
    ! initial temperature expressed in Kelvin
    iniTp = 2.0e4_num
    ! initial energy based on perfect gas law for fully ionised gas
    iniEn = ((2.0_num)/(gamma-1.0_num))*((kb_si*iniTp)/m0)
   

    ! Turn on or off the Braginskii thermal conduction term in the MHD equations   
    conduction = .FALSE.
    ! Apply a flux limiter to stop heat flows exceeding free streaming limit
    ! this is an experimental feature
    con_flux_limited = .FALSE.
    ! Fraction of free streaming heat flux used in limiter
    con_limit = 0.01_num
    ! Isotropic conductivity coefficient
    con_b_min = 0.001_num
    
    ! Use radiation as specified in SUBROUTINE calc_rad_losses in src/core/conrad.f90   
    radiation = .FALSE.
    ! TRUE if radiation losses are calculated implicitly along with conduction
    rad_implicit = .FALSE.
    ! Algorithm used to decide how radiation is calculated.
    ! 1. Rosner, R., Tucker, W. H., & Vaiana, G. S., 1978a, ApJ, 220, 643
    ! 2. Klimchuk, J. A., Patsourakos, S., & Cargill, P. J., 2008, ApJ, 1351
    ! First algorithm is used by default.
    rad_alg = 2
    ! The radiation loss threshold below which the timestep is not updated.
    rad_min_loss = 1.0e-12_num
    ! The change in temperature as a fraction of the actual temperature.
    rad_dTp_frac = 0.01_num
    ! True, if radiation is calculated for the average temperature over the timestep. 
    rad_avTp = .TRUE.

    ! Iteration count for Gauss-Seidel iteration with successive over-relaxation (SOR).
    ! This is the technique used to incorporate conduction and radiation.
    sor_it_lim = 1000
    ! Over-relaxation parameter.
    sor_w = 1.6_num
    ! Error threshold below which solution has converged.
    sor_err = 0.001_num
    ! Indicate how the absolute residual is weighted (default value is no weighting).
    ! 1: maximum energy (over entire grid) at time n, 2: current energy at time n, 3: iterated current energy.
    sor_err_wgt = 1

    ! True, if conduction and radiation are calculated over half the timestep. 
    conrad_half_dt = .TRUE.  

        
    ! loop radius
    lp_rad = 1.0_num

    ! the initial state for the alpha-omega configuration
    lp_alp = 1.8_num
    lp_omg = 3.0_num

    ! Turn on or off the initial (kink) perturbation
    perturbation = .TRUE.

    ! kink instability
    ptb_sign = 1.0_num
    ptb_k = 6.0_num*pi/20.0_num
    ptb_amp = 0.0001_num
    ptb_lim = 10.0_num
    
    ! Turn on or off lagrangian tracers, see lare3d.f90 and corks.f90
    diag_corks = .FALSE.
    
  END SUBROUTINE control_variables



  SUBROUTINE set_output_dumps

    ! The output directory for the code
    data_dir = "data"

    ! The interval between output snapshots. If SI_Input is true
    ! Then this is in seconds
    dt_snapshots = 5.0_num

    ! dump_mask is an array which specifies which quantities the
    ! code should output to disk in a data dump.
    ! The codes are
    ! 1  - rho
    ! 2  - energy
    ! 3  - vx
    ! 4  - vy
    ! 5  - vz
    ! 6  - bx
    ! 7  - by
    ! 8  - bz
    ! 9  - temperature
    ! 10 - pressure
    ! 11 - cs (sound speed)
    ! 12 - parallel_current
    ! 13 - perp_current
    ! 14 - neutral_faction
    ! 15 - eta_perp
    ! 16 - eta
    ! 17 - jx
    ! 18 - jy
    ! 19 - jz
    ! 20 - p_visc
    ! 21 - visc_heat
    ! 22 - visc_heat_xy
    ! 23 - visc_heat_xz
    ! 24 - visc_heat_yz
    ! 25 - visc_heat_xx
    ! 26 - visc_heat_yy
    ! 27 - visc_heat_zz
    ! 28 - ohmic_heat
    ! 29 - fx
    ! 30 - fy
    ! 31 - fz
    ! 32 - jcrit
    ! 33 - eta_num
    ! If a given element of dump_mask is true then that field is dumped
    ! If the element is false then the field isn't dumped
    ! N.B. if dump_mask(1:8) not true then the restart will not work
    dump_mask = .FALSE.
    dump_mask(1:10) = .TRUE. 
    dump_mask(17:33) = .TRUE. 
    IF (eos_number /= EOS_IDEAL)  dump_mask(14) = .TRUE. 
    IF (cowling_resistivity)  dump_mask(15) = .TRUE. 
    IF (resistive_mhd)  dump_mask(16) = .TRUE. 

  END SUBROUTINE set_output_dumps

END MODULE control
