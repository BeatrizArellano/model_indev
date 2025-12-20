module bio_params
  use precision_types,  only: rk
  use grids,            only: VerticalGrid
  use read_config_yaml, only: ConfigParams

  implicit none
  private

    public :: BioParams, read_bio_parameters
    public :: SedParams, read_sed_parameters
    public :: mol_diff

    type, public :: BioParams
        character(:), allocatable :: config_file          ! FABM configuration file
        logical  :: sediments_enabled = .false.
        logical  :: repair = .false.
        logical  :: output_conserved = .false.            ! Retrieve conserved quantities from FABM and output the column-integrated totals
        real(rk) :: frac_max 
        real(rk) :: cnpar 
        real(rk) :: min_dt
    end type BioParams

    type, public :: SedParams
        ! IMPORTANT (units):
        !   This data structure is used in two unit conventions, depending on where it is stored:
        !     - SedimentEnv%params_user : "literature/user units"
        !         * depths:       cm
        !         * burial rate:  cm yr-1
        !         * bioturbation: cm2 yr-1
        !         * irrigation:   yr-1
        !
        !     - SedimentEnv%params_SI   : "internal SI units"
        !         * depths:       m
        !         * burial rate:  m s-1
        !         * bioturbation: m2 s-1
        !         * irrigation:   s-1
        !
        !   Conversion is performed once during init_sediments() via convert_units_to_SI().
        !   Internal sediment computations must use params_SI only.
        !-----------------------------------------------------------------------------------
        !--- Burial rates
        real(rk) :: sed_rate               ! Linear sedimentation rate (user: cm/yr; SI: m/s)
        real(rk) :: sed_ref_depth          ! Depth in the sediments for the sedimentation rate measurement (user: cm;   SI: m)
        !--- Porosity
        real(rk) :: poro_sfc                ! Porosity at SWI [-]
        real(rk) :: poro_deep               ! Porosity at depth [-]
        real(rk) :: poro_decay              ! Decay depth for porosity exponential decay (user: cm; SI: m)
        !--- Bioturbation
        real(rk) :: biot_db_sfc             ! Bioturbation diffusivity at the sediment surface (user: cm2/yr; SI: m2/s)
        real(rk) :: biot_mld                ! Thickness of Sediment mixed layer depth below which bioturbation decreases exponentially (user: cm; SI: m)
        real(rk) :: biot_ez                 ! Coefficient (decay depth) for exponential bioturbation decrease (user: cm; SI: m)
        !--- Bioirrigation
        real(rk) :: irr_sfc                ! Exchange rate at the sediment-water interface (user: 1/yr; SI: 1/s) (alpha0 in Aller's model)
        real(rk) :: irr_ez                 ! Decay depth for exponential attenuation of irrigation (user: cm; SI: m)
    end type SedParams

    real(rk), parameter :: mol_diff   = 1.0e-6_rk             ! Molecular diffusivity in the water column [m2 s-1]

    ! ----------------- Default values for parameters ----------------------------------------------------------
    real(rk), parameter :: def_frac_max    = 0.1_rk      ! Maximum fractional change allowed per main timestep per tracer
    real(rk), parameter :: def_min_dt      = 10.0_rk     ! Minimum time-step in the inner loop
    real(rk), parameter :: def_cnpar       = 0.5_rk      ! Degree of Implicitness when solving diffusive mixing [0-1]
    logical,  parameter :: def_repair      = .false.     ! Indicates FABM whether to repair the state of the variables
    logical,  parameter :: def_conserv     = .false.     ! Retrieve conserved quantities from FABM and output the column-integrated totals
    logical,  parameter :: def_sed_enabled = .false.     ! Sediments enabled
    !-------------------------------------------------------------------------------------------------------------------

    ! ----------------- Default values for sediment parameters ----------------------------------------------------------
    ! --- Values are later converted to SI units (meters and seconds) where needed
    real(rk), parameter :: def_sed_rate   = 0.1_rk        ! Sedimentation rate. [cm yr-1]
    real(rk), parameter :: def_poro_sfc   = 0.95_rk       ! Porosity at the sediment–water interface. (Soetaert et al., 1996)
    real(rk), parameter :: def_poro_deep  = 0.80_rk       ! Porosity at depth in sediment. (Soetaert et al., 1996)
    real(rk), parameter :: def_poro_decay = 4.0_rk        ! Depth scale over which sediment porosity decreases exponentially with depth [cm] (Soetaert et al., 1996)
    real(rk), parameter :: def_db_sfc     = 5.0_rk        ! Bioturbation at the sediment surface [cm2 yr-1]
    real(rk), parameter :: def_biot_mld   = 5.0_rk        ! Sediment mixed layer depth below which bioturbation decreases exponentially [cm]
    real(rk), parameter :: def_biot_ez    = 1.0_rk        ! Coefficient (decay depth) for exponential bioturbation decrease [cm]
    real(rk), parameter :: def_irr_sfc    = 200_rk        ! Irrigation rate at the sediment-water interface [yr-1]
    real(rk), parameter :: def_irr_ez     = 2.0_rk        ! Decay depth for exponential attenuation of irrigation [cm]
    !-------------------------------------------------------------------------------------------------------------


contains
    ! Read YAML values (with validation) into a BioParams object.
    ! Any missing entries fall back to the module defaults.
    subroutine read_bio_parameters(cfg_params, bio)
        type(ConfigParams),  intent(in)   :: cfg_params
        type(BioParams), intent(out)      :: bio
        bio = default_bio_params()
        bio%config_file = cfg_params%get_param_str('biogeochemistry.config_file', default='fabm.yaml', trim_value=.true.)
        ! ---------------- Flags ----------------
        bio%sediments_enabled = cfg_params%get_param_logical('biogeochemistry.sediments.enabled', default=def_sed_enabled)
        bio%repair  = cfg_params%get_param_logical('biogeochemistry.repair_state', default=def_repair)
        bio%output_conserved  = cfg_params%get_param_logical('biogeochemistry.output_conserved_qt', default=def_conserv)
        ! ---------------- Mixing ----------------
        bio%cnpar = cfg_params%get_param_num('biogeochemistry.vertical_mixing.cnpar', default=def_cnpar, finite=.true., min=0._rk, max=1._rk)
        ! ---------------- Numerics ----------------
        bio%frac_max = cfg_params%get_param_num('biogeochemistry.numerics.max_change', default=def_frac_max, finite=.true., positive=.true.)
        bio%min_dt   = cfg_params%get_param_num('biogeochemistry.numerics.min_timestep', default=def_min_dt, finite=.true., positive=.true.)
    end subroutine read_bio_parameters


    ! Return the data structure pre-filled with defaults
    pure function default_bio_params() result(p)
        type(BioParams) :: p
        p%sediments_enabled = def_sed_enabled
        p%repair            = def_repair
        p%output_conserved  = def_conserv
        p%frac_max          = def_frac_max    
        p%cnpar             = def_cnpar
        p%min_dt            = def_min_dt
    end function default_bio_params


    ! Reading parameters relevant for sediments in units commonly reported in the literature
    subroutine read_sed_parameters(cfg_params, sed_grid, SedP)
        type(ConfigParams),       intent(in) :: cfg_params
        type(VerticalGrid),       intent(in) :: sed_grid          ! Sediment grid
        type(SedParams),          intent(inout) :: SedP

        real(rk) :: def_sed_ref_depth
        def_sed_ref_depth = sed_grid%depth * 100.0_rk  ! Default value for depth [cm] where sedimentation rate was measured
        ! --------- Sedimentation rate -----------
        SedP%sed_rate      = cfg_params%get_param_num('biogeochemistry.sediments.sedimentation_rate', default=def_sed_rate, finite=.true., min=0.0_rk)
        SedP%sed_ref_depth = cfg_params%get_param_num('biogeochemistry.sediments.sed_rate_depth', default=def_sed_ref_depth, finite=.true., min=0.0_rk)
        ! ---------------- Porosity profile ----------------
        SedP%poro_sfc   = cfg_params%get_param_num('biogeochemistry.sediments.porosity_surface', default=def_poro_sfc, finite=.true., min=0.0_rk, max=1._rk)
        SedP%poro_deep  = cfg_params%get_param_num('biogeochemistry.sediments.porosity_depth', default=def_poro_deep, finite=.true., min=0.0_rk, max=1._rk)
        SedP%poro_decay = cfg_params%get_param_num('biogeochemistry.sediments.porosity_decay_depth', default=def_poro_decay, finite=.true., min=0.0_rk)
        ! -------- Bioturbation profile ------------
        SedP%biot_db_sfc = cfg_params%get_param_num('biogeochemistry.sediments.db_surface', default=def_db_sfc, finite=.true., min=0.0_rk)
        SedP%biot_mld    = cfg_params%get_param_num('biogeochemistry.sediments.bioturbation_depth', default=def_biot_mld, finite=.true., min=0.0_rk)
        SedP%biot_ez     = cfg_params%get_param_num('biogeochemistry.sediments.bioturbation_decay_depth', default=def_biot_ez, finite=.true., positive=.true.)
        !--------- Bioirrigation profile -------
        SedP%irr_sfc = cfg_params%get_param_num('biogeochemistry.sediments.irrigation_surface', default=def_irr_sfc, finite=.true., min=0.0_rk)
        SedP%irr_ez  = cfg_params%get_param_num('biogeochemistry.sediments.irrigation_decay_depth', default=def_irr_ez, finite=.true., positive=.true.)
    end subroutine read_sed_parameters

end module bio_params