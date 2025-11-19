! src/physics/physics_main.F90
module physics_main
  use iso_fortran_env,     only: error_unit
  use precision_types,     only: rk, lk
  use read_config_yaml,    only: ConfigParams
  use geo_utils,           only: LocationInfo
  use grids,               only: VerticalGrid
  use tidal_readers,       only: read_tidal_parameters
  use tidal,               only: TidalSet, create_tidal_set, tide_pressure_slopes
  use forcing_manager,     only: ForcingSnapshot
  use EOS_eqns,            only: eos_density
  use momentum_eqns
  use turbulence,          only: TURBULENCE_ke, tke_min
  use heat_fluxes,         only: SURFACE_HEAT
  use vertical_mixing,     only: scalar_diffusion
  use numerical_stability, only: compute_phys_subcycles
  use tridiagonal,         only: init_tridiag, reset_tridiag
  use physics_params,      only: read_physics_parameters, &
                                 gravity, z0s_min, kappa, mol_nu, mol_diff_T
  use physics_types,       only: PhysicsState, PhysicsEnv
  use nan_checks,          only: check_nan_physics
  
  implicit none
  private

  ! Public API
  public :: init_physics, solve_physics, end_physics             ! set up grid, parameters, state

contains

    !======================
    ! Initialization
    !======================
    subroutine init_physics(cfg_params,location, grid, PE)
      ! Receives from main: user config and location
      type(ConfigParams), intent(in)    :: cfg_params
      type(LocationInfo), intent(in)    :: location
      type(VerticalGrid), intent(in)    :: grid
      type(PhysicsEnv),   intent(inout) :: PE

      integer :: k
      real(rk), allocatable :: height(:)
      real(rk) :: depth, zeta

  
      if (PE%is_init) return
      
      !---- Read Physics parameters and initial values
      call read_physics_parameters(cfg_params, PE%params)
  
      !--- Reading  tidal parameters ----
      call read_tidal_parameters(cfg_params, location%lat,location%lon, PE%Tides)   ! Reads tidal parameters
      call create_tidal_set(PE%Tides,location%lat)

      PE%grid = grid

      PE%PS%N = PE%grid%nz                    ! Number of layers in the water column, inside physics state    
      ! ---------- Allocating arrays and setting initial state --------------------------
      ! --- Arrays for prognostic variables (1..N) ---
      allocate(PE%PS%temp(PE%PS%N), PE%PS%sal(PE%PS%N), PE%PS%rho(PE%PS%N))
      allocate(PE%PS%velx(PE%PS%N), PE%PS%vely(PE%PS%N))
      PE%PS%temp = PE%params%temp0
      PE%PS%sal  = PE%params%sal0
      PE%PS%rho  = eos_density(PE%PS%temp, PE%PS%sal)
      PE%PS%velx = 0.0_rk
      PE%PS%vely = 0.0_rk
  
      ! --- Arrays for turbulence/mixing (N+1) ---
      allocate(PE%PS%Kz(0:PE%PS%N), PE%PS%Nz(0:PE%PS%N), PE%PS%tke(0:PE%PS%N), PE%PS%eps(0:PE%PS%N))
      allocate(PE%PS%cmue1(0:PE%PS%N))
      PE%PS%Kz  = 0.5_rk * PE%params%vismax     ! S2P3 style start
      PE%PS%Nz  = 0.5_rk * PE%params%vismax
      PE%PS%tke = tke_min
      PE%PS%eps = 5.0e-10_rk
      PE%PS%cmue1 = 0.0_rk
  
      ! --- A sensible initial length-scale (S2P3) ---
      depth = PE%grid%depth
      allocate(PE%PS%Lscale(0:PE%PS%N))
      allocate(height(0:PE%PS%N))
      ! Height above bed at interfaces from dz only
      height(0) = 0.0_rk
      do k = 1, PE%PS%N
        height(k) = height(k-1) + PE%grid%dz(k)
      end do

      ! Neutral mixing length 
      PE%PS%Lscale = 0.0_rk
      do k = 1, PE%PS%N-1
        zeta = height(k)/depth
        PE%PS%Lscale(k) = kappa * height(k) * sqrt(1.0_rk - zeta)
      end do
      ! Boundaries: copy neighbours
      PE%PS%Lscale(0)       = PE%PS%Lscale(1)
      PE%PS%Lscale(PE%PS%N) = PE%PS%Lscale(PE%PS%N-1)

      deallocate(height)

      ! Tridiagonal workspace
      call init_tridiag(PE%trid, PE%PS%N)
  
      PE%is_init = .true.
    end subroutine init_physics 

    !==================================================================================
    !                  Solves physics for each main time-step.
    !               It has an inner loop for numerical stability
    !==================================================================================
    subroutine solve_physics(PE, FS, dt_main, elapsed_time, is_first_step, ierr, errmsg)
      type(PhysicsEnv), intent(inout) :: PE
      type(ForcingSnapshot), intent(in)  :: FS
      integer(lk),        intent(in)  :: dt_main, elapsed_time
      logical,            intent(in)  :: is_first_step
      integer,           intent(out), optional :: ierr
      character(len=*),  intent(out), optional :: errmsg    

      integer  :: N, ti, n_sub
      real(rk) :: dtm, dt_sub, t_sub
      real(rk) :: Pxsum, Pysum  
      
      ! --- Local Arrays ---
      real(rk), allocatable :: u_old(:), u_new(:), v_old(:), v_new(:)
      real(rk), allocatable :: Kz_T(:), Nz_tot(:)

      dtm = real(dt_main, kind=rk)

      N = PE%PS%N                          ! Number of layers in the vertical grid
      if (allocated(u_old)) deallocate(u_old, u_new, v_old, v_new)
      allocate(u_old(N), u_new(N), v_old(N), v_new(N))
      if (allocated(Nz_tot)) deallocate(Nz_tot, Kz_T)   
      allocate(Nz_tot(0:N), Kz_T(0:N))  ! Arrays for viscosity and diffusivity plus molecular values

      ! Initialising them with current state
      u_old = PE%PS%velx;  v_old = PE%PS%vely

      
      ! Calculating wind stress from wind-speed and direction
      call wind_stress_from_speed_dir(FS%wind_spd, FS%wind_dir, PE%PS%tau_x, PE%PS%tau_y)
      ! Initial calculation for surface friction velocity and roughness length
      PE%PS%u_taus = (((PE%PS%tau_x/PE%PS%rho(N))**2 + (PE%PS%tau_y/PE%PS%rho(N))**2 ))**0.25_rk
      PE%PS%z0s = max(PE%params%charnock * (PE%PS%u_taus*PE%PS%u_taus) / gravity, z0s_min)  ! Surface roughness length

      ! Turbulence: once per main step     
        ! Solves turbulence using the Canuto k-ε closure scheme and calculates Kz and Nz
        call TURBULENCE_ke(N, dtm, PE%params, PE%grid%dz,                       &
                          density = PE%PS%rho, velx = u_old, vely = v_old,         &
                          u_taus = PE%PS%u_taus, u_taub = PE%PS%u_taub,            &
                          z0s = PE%PS%z0s, z0b = PE%PS%z0b,                        &
                          Kz = PE%PS%Kz, Nz = PE%PS%Nz, tke = PE%PS%tke,           &
                          eps = PE%PS%eps, Lscale=PE%PS%Lscale,                    &
                          cmue1=PE%PS%cmue1, trid=PE%trid, is_first_step=is_first_step)

        Nz_tot = PE%PS%Nz + mol_nu
        Kz_T   = PE%PS%Kz + mol_diff_T 

      ! Calculates the optimum size of the time-steps  for the inner loop
      ! by obtaining the stability condition so that explicit vertical diffusion of momentum is stable
      call compute_phys_subcycles(PE%grid%dz, Nz_tot, Kz_T, PE%params%vismax, dtm, n_sub, dt_sub, ierr, errmsg)    

      ! Inner subcycle to maintain numeriacl stability
      do ti = 1, n_sub
          t_sub = real(elapsed_time,rk) + (real(ti,rk) - 0.5_rk)*dt_sub   ! Elapsed time for each substep
  
          call SURFACE_HEAT(temp = PE%PS%temp, dz = PE%grid%dz, N = N, dt=dt_sub, &
                           rsds=FS%short_rad, rlds_down=FS%long_rad, wind_speed=FS%wind_spd,     &
                           airT=FS%air_temp, rh=FS%rel_hum, airP=FS%slp,                      &
                           lw_skin_penetration=PE%params%lw_skin_penetration,      &
                           lambda=PE%params%lambda)
          
          ! Update density
          PE%PS%rho  = eos_density(PE%PS%temp, PE%PS%sal)
          

  
          ! Calculate current profile from tidal currents
          call tide_pressure_slopes(PE%Tides, t_sub, Pxsum, Pysum)             ! Pressure-gradient slopes from tidal constituents at this time
  
          ! Accelerate by pressure gradients (x and y components)
          call EQN_PRESSURE(dt_sub, Pxsum, u_old)  ! Already updates u_old inplace
          call EQN_PRESSURE(dt_sub, Pysum, v_old)


  
          ! Apply surface/bottom stresses and vertical viscosity (x component)
          call EQN_FRICTION( vel_comp_old = u_old, vel_comp_new = u_new,                &
                             vel_comp2_bottom = v_old(1), Nz=Nz_tot, h=PE%grid%dz,      &
                             dt=dt_sub, h0b=PE%params%h0b, density=PE%PS%rho,           &
                             tau_surf=PE%PS%tau_x, tau_surf2=PE%PS%tau_y,               &
                             charnock=PE%params%charnock,                               &
                             u_taub=PE%PS%u_taub, z0b=PE%PS%z0b, stressb=PE%PS%stressb, &
                             u_taus=PE%PS%u_taus, z0s=PE%PS%z0s )
  
          ! Apply surface/bottom stresses and vertical viscosity (y component)
          call EQN_FRICTION( vel_comp_old = v_old, vel_comp_new = v_new,                &
                             vel_comp2_bottom = u_old(1), Nz=Nz_tot, h=PE%grid%dz,      &
                             dt=dt_sub, h0b=PE%params%h0b, density=PE%PS%rho,           &
                             tau_surf=PE%PS%tau_y, tau_surf2=PE%PS%tau_x,               &
                             charnock=PE%params%charnock,                               &
                             u_taub=PE%PS%u_taub, z0b=PE%PS%z0b, stressb=PE%PS%stressb, &
                             u_taus=PE%PS%u_taus, z0s=PE%PS%z0s )
          ! Update olds for next substep
          u_old = u_new;  v_old = v_new
  
          ! Rotate velocities due to the Coriolis force
          call EQN_CORIOLIS(u_old, v_old, dt_sub, PE%Tides%f0)
  
          ! Update the state of velocities
          PE%PS%velx = u_old
          PE%PS%vely = v_old
  
          ! Then mix the vertical thermal structure with semi-implicit scalar scheme      
          call scalar_diffusion(PE%PS%temp, N, dt_sub, PE%grid%dz, Kz_T, PE%params%cnpar, PE%trid, ierr)
          ! Recompute density after thermal mixing
          PE%PS%rho = eos_density(PE%PS%temp, PE%PS%sal)

          call check_nan_physics(PE%PS)
      end do
      
!write(*,*) 'u_taub=', PE%PS%u_taub, ' Kz(B)=', PE%PS%Kz(0), ' tke=', PE%PS%tke(0)
      deallocate(u_old, u_new, v_old, v_new)
    end subroutine solve_physics  
  
    !==============================
    ! Clear Physics Environment
    !===============================
    subroutine end_physics(PE)    
      type(PhysicsEnv),   intent(inout) :: PE

      if (allocated(PE%PS%temp))   deallocate(PE%PS%temp, PE%PS%sal, PE%PS%rho)
      if (allocated(PE%PS%velx))   deallocate(PE%PS%velx, PE%PS%vely)
      if (allocated(PE%PS%Kz))     deallocate(PE%PS%Kz, PE%PS%Nz, PE%PS%tke, PE%PS%eps)
      if (allocated(PE%PS%Lscale)) deallocate(PE%PS%Lscale)
      if (allocated(PE%PS%cmue1))  deallocate(PE%PS%cmue1)
      PE%PS%N = 0
      ! reset tides (deallocates allocatable components via intrinsic assignment)
      PE%Tides = TidalSet()
      call reset_tridiag(PE%trid)    

      PE%is_init = .false.        
   end subroutine end_physics

end module physics_main
