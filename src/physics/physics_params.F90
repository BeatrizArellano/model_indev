module physics_params
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams

  implicit none
  private

  public :: gravity, rho0, rho_air, cp_sw, cp_air, kappa, sigma_SB, Omega, p_atm_dbar
  public :: mol_vis, mol_diff_T, mol_diff_S, mol_nu, z0s_min

  ! ------------------ Constants ----------------------------------
  ! ---- Universal constants ----
  real(rk), parameter :: gravity         = 9.81_rk               ! Gravity acceleration
  real(rk), parameter :: Omega           = 7.2921159e-5_rk       ! Earth's angular velocity [rad/s]
  real(rk), parameter :: rho0            = 1025.0_rk             ! Reference seawater density [kg m-3]
  real(rk), parameter :: rho_air         = 1.3_rk                ! Reference air density [kg m-3] Using 1.3 as in S2P3
  real(rk), parameter :: cp_sw           = 3900.0_rk             ! seawater heat capacity [J kg-1 K-1]
  real(rk), parameter :: cp_air          = 1004.0_rk             ! air heat capacity [J kg-1 K-1]
  real(rk), parameter :: kappa           = 0.41_rk               ! von Kármán constant
  real(rk), parameter :: sigma_SB        = 5.670374419e-8_rk     ! Stefan-Boltzman 
  real(rk), parameter :: p_atm_dbar      = 10.1325_rk            ! Atmospheric pressure at sea level 1 atm [dbar]
  !----- Specific constants for viscosity and diffusivity
  real(rk), parameter :: mol_vis         = 1.0e-6_rk             ! Molecular viscosity [m2 s-1]
  real(rk), parameter :: mol_diff_T      = 1.4e-7_rk             ! Thermal molecular diffusivity [m2 s-1]
  real(rk), parameter :: mol_diff_S      = 1.0e-9_rk             ! Salt molecular diffusivity [m2 s-1]  
  real(rk), parameter :: mol_nu          = 1.3e-6_rk             ! Molecular (kinematic) viscosity of water [m²/s]
  real(rk), parameter :: z0s_min         = 0.02_rk               ! Minimum value for hydrodynamic roughness [m]


  ! ----------------- Default values for parameters ----------------------------------------------------------
  real(rk), parameter :: def_temp0               = 12.0_rk       ! Initial temperature
  real(rk), parameter :: def_sal0                = 35.0_rk       ! Constant salinity
  real(rk), parameter :: def_charnock            = 1400.0_rk     ! Empirical constant for roughness adaptation (Charnock, 1955)
  real(rk), parameter :: def_h0b                 = 0.05_rk       ! Bottom roughness height (h0b)
  real(rk), parameter :: def_vismax              = 0.1_rk        ! Maximum diffusivity and viscosity
  real(rk), parameter :: def_Kz_bg               = 1.0e-5_rk     ! Background diffusivity [m2 s-1]
  real(rk), parameter :: def_Nz_bg               = 1.0e-5_rk     ! Background viscosity [m2 s-1]
  real(rk), parameter :: def_frac_nonvis         = 0.58_rk       ! Fraction of SWR in the non-visible band (UV+IR)
  real(rk), parameter :: def_depth_nonvis        = 0.35_rk       ! e-folding depth of non-visible band [m]
  real(rk), parameter :: def_depth_vis           = 20.0_rk       ! e-folding depth of visible band [m]
  real(rk), parameter :: def_frac_par            = 0.45_rk       ! Fraction of shortwave radiation treated as Photosynthetically Active Radiation (PAR)
  real(rk), parameter :: def_lw_skin_penetration = 0.95_rk       ! Fraction of incoming LW that reaches below the skin 
  real(rk), parameter :: def_cnpar               = 0.5_rk        ! Degree of Implicitness when solving diffusive mixing [0-1]
  character(len=8), parameter :: def_salmode = 'constant'        ! Default mode for salinity
  !--------------------------------------------------------------------------------------------------------------------

  public :: read_physics_parameters
  public :: PhysicsParams

  ! Container for all physics parameters
  type, public :: PhysicsParams
     ! Initial conditions
     real(rk) :: temp0               ! [°C]
     real(rk) :: sal0                ! [PSU]
     ! Compute salinity
     logical  :: compute_salinity  
     ! Surface
     real(rk) :: charnock
     ! Seabed         
     real(rk) :: h0b                 ! roughness height [m]
     ! Mixing         
     real(rk) :: vismax              ! [m2 s-1]
     real(rk) :: Kz_bg               ! [m2 s-1]
     real(rk) :: Nz_bg               ! [m2 s-1]
     ! Heat and radiation
     real(rk) :: frac_nonvis         ! [-]
     real(rk) :: depth_nonvis        ! [m]
     real(rk) :: depth_vis           ! [m]
     real(rk) :: frac_par        ! [-]
     real(rk) :: lw_skin_penetration ! [fraction 0-1]
     logical  :: apply_heat_bioshade
     logical  :: apply_par_bioshade

     ! Numerics
     real(rk) :: cnpar               ! implicitness [0..1]
  end type PhysicsParams

contains

    ! Read YAML values (with validation) into a PhysicsParams object.
    ! Any missing entries fall back to the module defaults.
    subroutine read_physics_parameters(cfg_params, phys)
        type(ConfigParams),  intent(in)   :: cfg_params
        type(PhysicsParams), intent(out)  :: phys

        logical :: is_bio_enabled
        character(len=8), dimension(2) :: sal_choices
        character(:), allocatable :: sal_mode
        sal_choices = ['constant','compute ']

        

        phys = default_physics_params()
        ! ---------------- Variables ----------------
        phys%temp0 = cfg_params%get_param_num('physics.variables.temperature.initial_value', default=def_temp0, finite=.true.)
        phys%sal0  = cfg_params%get_param_num('physics.variables.salinity.initial_value', default=def_sal0, finite=.true., min=0._rk)
        ! -------------- Flag
        sal_mode  = cfg_params%get_param_str('physics.variables.salinity.mode', &
                                             default=def_salmode, choices=sal_choices, trim_value=.true., match_case=.false.)
        if (sal_mode == 'compute') then
            phys%compute_salinity = .true.
        else 
            phys%compute_salinity = .false.
        end if
        !------------------Surface ----------------
        phys%charnock = cfg_params%get_param_num('physics.surface.roughness.charnock', default=def_charnock, finite=.true., min=0._rk)
        ! ---------------- Seabed ----------------
        phys%h0b   = cfg_params%get_param_num('physics.seabed.roughness_height', default=def_h0b, finite=.true., min=0._rk)
        ! ---------------- Mixing ----------------
        phys%vismax = cfg_params%get_param_num('physics.mixing.vismax', default=def_vismax, finite=.true., min=0._rk)
        phys%Kz_bg  = cfg_params%get_param_num('physics.mixing.Kz_bg',  default=def_Kz_bg,  finite=.true., min=0._rk)
        phys%Nz_bg  = cfg_params%get_param_num('physics.mixing.Nz_bg',  default=def_Nz_bg,  finite=.true., min=0._rk)
        ! ---------------- Heat and radiation ----------------
        phys%frac_nonvis        = cfg_params%get_param_num('physics.radiation.nonvisible_fraction', default=def_frac_nonvis, finite=.true., min=0._rk, max=1.0_rk)
        phys%depth_nonvis       = cfg_params%get_param_num('physics.radiation.depth_nonvisible', default=def_depth_nonvis, finite=.true., positive=.true.)       
        phys%depth_vis          = cfg_params%get_param_num('physics.radiation.depth_visible', default=def_depth_vis, finite=.true., positive=.true.)   
        phys%frac_par           = cfg_params%get_param_num('physics.radiation.par_fraction', default=def_frac_par, finite=.true., min=0._rk, max=1.0_rk)
        phys%lw_skin_penetration = cfg_params%get_param_num('physics.radiation.lw_skin_penetration', default=def_lw_skin_penetration, finite=.true., min=0._rk, max=1._rk)

        phys%apply_heat_bioshade = cfg_params%get_param_logical('physics.radiation.apply_bioshading_to_heat', default=.true.)
        phys%apply_par_bioshade  = cfg_params%get_param_logical('physics.radiation.apply_bioshading', default=.true.)

        is_bio_enabled = cfg_params%get_param_logical('biogeochemistry.enabled', default=.false.)
        if (.not. is_bio_enabled) then
            if(phys%apply_heat_bioshade) write(*,*) 'Warning: apply_bioshading_to_heat set to "no" since Biogeochemistry is disabled.'
            phys%apply_heat_bioshade = .false. 
            phys%apply_par_bioshade  = .false.
        end if
        ! ---------------- Numerics ----------------
        phys%cnpar = cfg_params%get_param_num('physics.mixing.cnpar', default=def_cnpar, finite=.true., min=0._rk, max=1._rk)
    end subroutine read_physics_parameters


    ! Return a struct pre-filled with defaults
    pure function default_physics_params() result(p)
        type(PhysicsParams) :: p
        p%charnock = def_charnock
        p%compute_salinity = .false.   ! consistent with def_salmode='constant'
        p%temp0      = def_temp0
        p%sal0       = def_sal0
        p%h0b        = def_h0b
        p%vismax     = def_vismax
        p%Kz_bg      = def_Kz_bg
        p%Nz_bg      = def_Nz_bg
        p%frac_nonvis = def_frac_nonvis
        p%depth_nonvis = def_depth_nonvis
        p%depth_vis   = def_depth_vis
        p%frac_par   = def_frac_par
        p%lw_skin_penetration = def_lw_skin_penetration
        p%apply_heat_bioshade = .true.
        p%apply_par_bioshade  = .true.
        p%cnpar      = def_cnpar
    end function default_physics_params


end module physics_params