module bio_params
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams

  implicit none
  private

  public :: BioParams, read_bio_parameters
  public :: mol_diff

  real(rk), parameter :: mol_diff   = 1.0e-6_rk             ! Molecular diffusivity in the water column [m2 s-1]

  ! ----------------- Default values for parameters ----------------------------------------------------------
  real(rk), parameter :: def_frac_max    = 0.1_rk      ! Maximum fractional change allowed per main timestep per tracer
  real(rk), parameter :: def_min_dt      = 10.0_rk     ! Minimum time-step in the inner loop
  real(rk), parameter :: def_cnpar       = 0.5_rk      ! Degree of Implicitness when solving diffusive mixing [0-1]
  logical,  parameter :: def_repair      = .false.     ! Indicates FABM whether to repair the state of the variables
  logical,  parameter :: def_conserv     = .false.     ! Retrieve conserved quantities from FABM and output the column-integrated totals
  logical,  parameter :: def_sed_enabled = .false.     ! Sediments enabled
  !--------------------------------------------------------------------------------------------------------------------

type, public :: BioParams
    character(:), allocatable :: config_file          ! FABM configuration file
    logical  :: sediments_enabled = .false.
    logical  :: repair = .false.
    logical  :: output_conserved = .false.            ! Retrieve conserved quantities from FABM and output the column-integrated totals
    real(rk) :: frac_max 
    real(rk) :: cnpar 
    real(rk) :: min_dt
end type BioParams


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


    ! Return a struct pre-filled with defaults
    pure function default_bio_params() result(p)
        type(BioParams) :: p
        p%sediments_enabled = def_sed_enabled
        p%repair            = def_repair
        p%output_conserved  = def_conserv
        p%frac_max          = def_frac_max    
        p%cnpar             = def_cnpar
        p%min_dt            = def_min_dt
    end function default_bio_params

end module bio_params