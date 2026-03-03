! src/physics/physics_main.F90
!=======================================================================================
!   Shelf Sea Physics:
!   1-D MODEL OF THE EQUATION OF MOTION USING THE Canuto k-e TURBULENCE CLOSURE SCHEME
!
!   This code is fully based on the S2P3 model lineage:
!
!   • Original S2P3 (v7.0):
!       Jonathan Sharples – Univ. of Liverpool & NOC
!       Documented in: Simpson & Sharples (2012), CUP.
!
!   • Regional framework S2P3-R (v1.0):
!       Marsh, Hickman & Sharples (2015), GMD 8, 3163–3178.
!
!   • Large-scale/efficient S2P3-R v2.0:
!       Halloran et al. (2021), GMD 14, 6177–6195.
!
!=======================================================================================
module physics_main  
  use EOS_eqns,            only: eos_density
  use forcing_manager,     only: ForcingSnapshot
  use freshwater_fluxes,   only: apply_surface_freshwater
  use geo_utils,           only: LocationInfo
  use grids,               only: VerticalGrid
  use heat_fluxes,         only: compute_heat_tendency, apply_temperature_tendency
  use momentum_eqns,       only: EQN_PRESSURE, EQN_CORIOLIS,   &
                                 compute_surface_stress, compute_bottom_stress
  use nan_checks,          only: check_nan_physics                               
  use numerical_stability, only: compute_phys_substeps
  use physics_params,      only: read_physics_parameters, &
                                 gravity, kappa, mol_nu, mol_diff_T, rho0
  use physics_types,       only: PhysicsState, PhysicsEnv       
  use phys_var_registry,   only: register_physics_variables    
  use precision_types,     only: rk, lk               
  use radiation,           only: compute_par_profile    
  use read_config_yaml,    only: ConfigParams 
  use tidal_readers,       only: read_tidal_parameters
  use tidal,               only: TidalSet, create_tidal_set, tide_pressure_accel
  use tridiagonal,         only: init_tridiag, clear_tridiag
  use turbulence,          only: TURBULENCE_ke, tke_min
  use vertical_mixing,     only: scalar_diffusion, BC_NEUMANN  
  
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
      type(VerticalGrid), intent(in)    :: grid     ! water grid
      type(PhysicsEnv),   intent(inout) :: PE

      integer :: k, N
      real(rk), allocatable :: height(:)
      real(rk) :: depth, zeta

      if (PE%is_init) return
      
      !---- Read Physics parameters and initial values
      call read_physics_parameters(cfg_params, PE%params)
  
      !--- Reading  tidal parameters ----
      call read_tidal_parameters(cfg_params, location%lat,location%lon, PE%Tides)   ! Reads tidal parameters
      call create_tidal_set(PE%Tides,location%lat)

      PE%grid = grid

      N       = PE%grid%nz
      PE%PS%N = PE%grid%nz                    ! Number of layers in the water column, inside physics state    
      ! ---------- Allocating arrays and setting initial state --------------------------
      ! --- Arrays for prognostic variables (1..N) ---
      allocate(PE%PS%temp(PE%PS%N), PE%PS%sal(PE%PS%N), PE%PS%rho(PE%PS%N))
      allocate(PE%PS%velx(PE%PS%N), PE%PS%vely(PE%PS%N))
      allocate(PE%PS%swr(PE%PS%N), PE%PS%par(PE%PS%N))
      PE%PS%temp = PE%params%temp0
      PE%PS%sal  = PE%params%sal0
      PE%PS%rho  = eos_density(PE%PS%temp, PE%PS%sal)
      PE%PS%velx = 0.0_rk
      PE%PS%vely = 0.0_rk
      PE%PS%swr  = 0.0_rk 
      PE%PS%par  = 0.0_rk 
  
      ! --- Arrays for turbulence/mixing (N+1) ---
      allocate(PE%PS%Kz(0:PE%PS%N), PE%PS%Nz(0:PE%PS%N), PE%PS%tke(0:PE%PS%N), PE%PS%eps(0:PE%PS%N))
      allocate(PE%PS%NN(0:PE%PS%N), PE%PS%SS(0:PE%PS%N), PE%PS%Ri(0:PE%PS%N))
      allocate(PE%PS%cmue1(0:PE%PS%N))
      PE%PS%Kz  = 0.5_rk * PE%params%vismax     ! Initial values like in S2P3
      PE%PS%Nz  = 0.5_rk * PE%params%vismax
      PE%PS%tke = tke_min
      PE%PS%eps = 5.0e-10_rk
      PE%PS%cmue1 = 0.0_rk

      !--- Allocate working arrays
      if (allocated(PE%u_old)) deallocate(PE%u_old, PE%u_new, PE%v_old, PE%v_new)
      allocate(PE%u_old(N), PE%u_new(N), PE%v_old(N), PE%v_new(N))
      if (allocated(PE%Nz_tot)) deallocate(PE%Nz_tot, PE%Kz_T)   
      allocate(PE%Nz_tot(0:N), PE%Kz_T(0:N))  ! Arrays for viscosity and diffusivity plus molecular values
      if (allocated(PE%dTdt_heat)) deallocate(PE%dTdt_heat)   
      allocate(PE%dTdt_heat(N))               ! Array for temperature tendencies computed from heat fluxes
      if (allocated(PE%atten_bio)) deallocate(PE%atten_bio)   
      allocate(PE%atten_bio(N))               ! Array for attenuation coefficients from biogeochemsitry
      PE%atten_bio = 0.0_rk      
  
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

      ! Allocate arrays for the tridiagonal solver
      call init_tridiag(PE%trid, PE%PS%N)

      ! Register variables and define the ones to include in the output
      call register_physics_variables(cfg_params, PE)
  
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

        ! Locals
        integer  :: N, ti, n_sub
        real(rk) :: dtm, dt_sub, t_sub
        real(rk) :: Pxsum, Pysum  
        real(rk) :: flux_u_bot, flux_v_bot
        real(rk) :: flux_u_sfc, flux_v_sfc

        integer :: ierr_l
        character(len=256) :: errmsg_l

        ierr_l = 0
        errmsg_l = ''

        dtm = real(dt_main, kind=rk)

        N = PE%PS%N                          ! Number of layers in the vertical grid

        ! Initialising them with current state
        PE%u_old = PE%PS%velx;  PE%v_old = PE%PS%vely

        ! Calculating wind stress, surface friction velocity and roughness length
        call compute_surface_stress(u10= FS%wind_u10, v10= FS%wind_v10, charnock= PE%params%charnock,   &
                                     tau_x= PE%PS%tau_x, tau_y= PE%PS%tau_y, u_taus= PE%PS%u_taus,      &
                                     z0s= PE%PS%z0s, windspeed = PE%PS%wind_speed)   
    
        ! Update density
        PE%PS%rho  = eos_density(PE%PS%temp, PE%PS%sal)

        call compute_heat_tendency(temp=PE%PS%temp, dz=PE%grid%dz, N=N,                                               &
                                   rsds=FS%short_rad, rlds_down=FS%long_rad, wind_speed=PE%PS%wind_speed,             &
                                   airT=FS%air_temp, rh=FS%rel_hum, airP=FS%slp,                                      &
                                   lw_skin_penetration=PE%params%lw_skin_penetration,                                 &
                                   nonvisible_fraction=PE%params%frac_nonvis, depth_nonvisible=PE%params%depth_nonvis, &
                                   depth_visible=PE%params%depth_vis, deposit_bottom_residual=.true.,                 &
                                   apply_bioshading_to_heat=PE%params%apply_heat_bioshade,                            &
                                   dTdt_heat=PE%dTdt_heat, Q_net_surf=PE%PS%Q_net_surf, max_abs_dTdt=PE%max_abs_dTdt, &
                                   swr_c=PE%PS%swr, atten_bio=PE%atten_bio)
        
        call compute_par_profile(N=N, dz=PE%grid%dz, rsds=FS%short_rad,                            &
                                 par_fraction=PE%params%frac_par, depth_visible=PE%params%depth_vis,                                 &
                                 apply_bioshading_to_par=PE%params%apply_par_bioshade, atten_bio=PE%atten_bio,                          &
                                 par=PE%PS%par, par_sfc=PE%PS%par_sfc)
 
        ! Turbulence: once per main step     
        ! Solves turbulence using the Canuto k-eps closure scheme and calculates Kz and Nz
        call TURBULENCE_ke(N, dtm, PE%params, PE%grid%dz,                                &
                          density = PE%PS%rho, velx = PE%u_old, vely = PE%v_old,         &
                          u_taus = PE%PS%u_taus, u_taub = PE%PS%u_taub,                  &
                          z0s = PE%PS%z0s, z0b = PE%PS%z0b,                              &
                          Kz = PE%PS%Kz, Nz = PE%PS%Nz, tke = PE%PS%tke,                 &
                          eps = PE%PS%eps, Lscale=PE%PS%Lscale,                          &
                          NN=PE%PS%NN, SS=PE%PS%SS, Ri=PE%PS%Ri,                         &
                          cmue1=PE%PS%cmue1, trid=PE%trid, is_first_step=is_first_step)

        PE%Nz_tot = PE%PS%Nz + mol_nu
        PE%Kz_T   = PE%PS%Kz + mol_diff_T 

        ! Compute bed stress components (N/m^2) 
        call compute_bottom_stress(u=PE%u_old, v=PE%v_old, dz_btm=PE%grid%dz(1), h0b=PE%params%h0b, &
                                  tau_bx=PE%PS%tau_bx, tau_by=PE%PS%tau_by, u_taub=PE%PS%u_taub, z0b=PE%PS%z0b, stressb=PE%PS%stressb)

        ! Calculates the optimum size of the time-steps  for the inner loop
        ! by obtaining the stability condition so that diffusion of properties and heat tendencies are stable
        call compute_phys_substeps(dz=PE%grid%dz, Nz=PE%Nz_tot, Kz=PE%Kz_T,                          &
                                  cnpar=PE%params%cnpar, dt_main=dtm,  max_abs_dTdt=PE%max_abs_dTdt, &
                                  tau_bx=PE%PS%tau_bx, tau_by=PE%PS%tau_by, tau_x=PE%PS%tau_x,       &        
                                  tau_y=PE%PS%tau_y, n_sub=n_sub,dt_sub=dt_sub,                      &                          
                                  ierr=ierr_l, errmsg= errmsg_l)

        ! Inner subcycle to maintain numeriacl stability
        do ti = 1, n_sub
            t_sub = real(elapsed_time,rk) + (real(ti,rk) - 0.5_rk)*dt_sub   ! Elapsed time for each substep   

            call apply_temperature_tendency(PE%PS%temp, N, dt_sub, PE%dTdt_heat)

            ! Then mix the vertical thermal structure with semi-implicit scalar scheme      
            call scalar_diffusion(PE%PS%temp, N, dt_sub, PE%grid%dz, PE%Kz_T, PE%params%cnpar, PE%trid, ierr_l)

            ! Recompute density after thermal mixing
            PE%PS%rho = eos_density(PE%PS%temp, PE%PS%sal)

            if (PE%params%compute_salinity) then
              ! Doing nothing at the moment
              call apply_surface_freshwater(salt=PE%PS%sal, dz=PE%grid%dz, precip=FS%precip, evap=FS%evap, runoff=FS%runoff)
              ! Update density
              !PE%PS%rho  = eos_density(PE%PS%temp, PE%PS%sal)
              ! Then mix freshwater in the water column
              !call scalar_diffusion(PE%PS%salt, N, dt_sub, PE%grid%dz, PE%Kz_T, PE%params%cnpar, PE%trid, ierr)
            end if            

    
            ! Calculate current profile from tidal currents
            call tide_pressure_accel(PE%Tides, t_sub, Pxsum, Pysum)             ! Pressure-gradient accelerations from tidal constituents at this time
    
            ! Accelerate by pressure gradients (x and y components)
            call EQN_PRESSURE(dt_sub, Pxsum, Pysum, PE%u_old, PE%v_old) ! Already updates u_old inplace

            ! Compute bed stress components (N/m^2) 
            call compute_bottom_stress(u=PE%u_old, v=PE%v_old, dz_btm=PE%grid%dz(1), h0b=PE%params%h0b, &
                                       tau_bx=PE%PS%tau_bx, tau_by=PE%PS%tau_by, u_taub=PE%PS%u_taub, z0b=PE%PS%z0b, stressb=PE%PS%stressb)

            ! Convert dynamic stresses (Pa) to kinematic momentum fluxes (m^2 s^-2)
            ! to use as Neumann boundary conditions in the momentum diffusion scheme
            flux_u_bot = PE%PS%tau_bx / rho0; flux_v_bot = PE%PS%tau_by / rho0
            flux_u_sfc  = PE%PS%tau_x  / rho0; flux_v_sfc  = PE%PS%tau_y  / rho0

            ! Apply surface/bottom stresses and vertical diffusivity of momentum (using a semi-implicit scheme)
            ! u_old and v_old are updated inplace
            call scalar_diffusion(PE%u_old, N, dt_sub, PE%grid%dz, PE%Nz_tot, PE%params%cnpar, PE%trid, ierr_l, &
                                  bc_bot_type=BC_NEUMANN, bc_bot_value=flux_u_bot, &
                                  bc_top_type=BC_NEUMANN, bc_top_value=flux_u_sfc)

            call scalar_diffusion(PE%v_old, N, dt_sub, PE%grid%dz, PE%Nz_tot, PE%params%cnpar, PE%trid, ierr_l, &
                                  bc_bot_type=BC_NEUMANN, bc_bot_value=flux_v_bot, &
                                  bc_top_type=BC_NEUMANN, bc_top_value=flux_v_sfc)

            PE%u_new = PE%u_old
            PE%v_new = PE%v_old
    
            ! Rotate velocities due to the Coriolis force
            call EQN_CORIOLIS(PE%u_old, PE%v_old, dt_sub, PE%Tides%f0)
    
            ! Update the state of velocities
            PE%PS%velx = PE%u_old
            PE%PS%vely = PE%v_old  

            call check_nan_physics(PE%PS)
        end do

      if (present(ierr))   ierr = ierr_l
      if (present(errmsg)) errmsg = trim(errmsg_l)

    end subroutine solve_physics  
  
    !==============================
    ! Clear Physics Environment
    !===============================
    subroutine end_physics(PE)    
      type(PhysicsEnv),   intent(inout) :: PE

      integer :: i
      if (allocated(PE%PS%temp))   deallocate(PE%PS%temp)
      if (allocated(PE%PS%sal))    deallocate(PE%PS%sal)
      if (allocated(PE%PS%rho))    deallocate(PE%PS%rho)

      if (allocated(PE%PS%swr))    deallocate(PE%PS%swr)
      if (allocated(PE%PS%par))    deallocate(PE%PS%par)

      if (allocated(PE%PS%velx))   deallocate(PE%PS%velx)
      if (allocated(PE%PS%vely))   deallocate(PE%PS%vely)

      if (allocated(PE%PS%Kz))     deallocate(PE%PS%Kz)
      if (allocated(PE%PS%Nz))     deallocate(PE%PS%Nz)
      if (allocated(PE%PS%tke))    deallocate(PE%PS%tke)
      if (allocated(PE%PS%eps))    deallocate(PE%PS%eps)
      if (allocated(PE%PS%NN))     deallocate(PE%PS%NN)   
      if (allocated(PE%PS%SS))     deallocate(PE%PS%SS)   
      if (allocated(PE%PS%Ri))     deallocate(PE%PS%Ri)   

      if (allocated(PE%PS%Lscale)) deallocate(PE%PS%Lscale)
      if (allocated(PE%PS%cmue1))  deallocate(PE%PS%cmue1)

      if (allocated(PE%u_old)) deallocate(PE%u_old)
      if (allocated(PE%u_new)) deallocate(PE%u_new)
      if (allocated(PE%v_old)) deallocate(PE%v_old)
      if (allocated(PE%v_new)) deallocate(PE%v_new)

      if (allocated(PE%Nz_tot))    deallocate(PE%Nz_tot)
      if (allocated(PE%Kz_T))      deallocate(PE%Kz_T)
      if (allocated(PE%dTdt_heat)) deallocate(PE%dTdt_heat)
      if (allocated(PE%atten_bio)) deallocate(PE%atten_bio)
      

      PE%PS%N = 0

      if (allocated(PE%phys_vars)) then
        do i = 1, size(PE%phys_vars)
          nullify(PE%phys_vars(i)%data_0d)
          nullify(PE%phys_vars(i)%data_1d)
        end do
        deallocate(PE%phys_vars)
      end if

      PE%Tides = TidalSet()
      call clear_tridiag(PE%trid)
   

      PE%is_init = .false.        
   end subroutine end_physics

end module physics_main
