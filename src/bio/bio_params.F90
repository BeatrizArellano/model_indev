module bio_params
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams

  implicit none
  private

  public :: BioParams, is_bio_enabled, read_bio_parameters
  public :: mol_diff

  real(rk), parameter :: mol_diff   = 1.0e-6_rk             ! Molecular diffusivity in the water column [m2 s-1]

  ! ----------------- Default values for parameters ----------------------------------------------------------
  real(rk), parameter :: def_frac_max    = 0.1_rk      ! Initial temperature
  real(rk), parameter :: def_cnpar       = 0.5         ! Degree of Implicitness when solving diffusive mixing [0-1]
  logical,  parameter :: def_repair      = .false.     ! Indicates FABM whether to repair the state of the variables
  logical,  parameter :: def_sed_enabled = .false.     ! Sediments enabled
  !--------------------------------------------------------------------------------------------------------------------

type, public :: BioParams
    character(:), allocatable :: config_file          ! FABM configuration file
    logical  :: sediments_enabled = .false.
    logical  :: repair = .false.
    real(rk) :: frac_max 
    real(rk) :: cnpar 
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
        ! ---------------- Mixing ----------------
        bio%cnpar = cfg_params%get_param_num('biogeochemistry.vertical_mixing.cnpar', default=def_cnpar, finite=.true., min=0._rk, max=1._rk)
        ! ---------------- Numerics ----------------
        bio%frac_max = cfg_params%get_param_num('biogeochemistry.numerics.max_change', default=def_frac_max, finite=.true., positive=.true.)
    end subroutine read_bio_parameters


    ! Return a struct pre-filled with defaults
    pure function default_bio_params() result(p)
        type(BioParams) :: p
        p%sediments_enabled = def_sed_enabled
        p%repair            = def_repair
        p%frac_max          = def_frac_max    
        p%cnpar             = def_cnpar
    end function default_bio_params
    
    subroutine is_bio_enabled(cfg_params, bio_enabled, sediments_enabled)
        type(ConfigParams),  intent(in)   :: cfg_params
        logical,             intent(out)  :: bio_enabled, sediments_enabled

        bio_enabled       = cfg_params%get_param_logical('biogeochemistry.enabled', default=.false.)
        sediments_enabled = cfg_params%get_param_logical('biogeochemistry.sediments.enabled', default=.false.)

        ! Enforce hierarchy: no sediments without column biogeochemistry
        if (.not. bio_enabled) then
            sediments_enabled = .false.
        end if
    end subroutine is_bio_enabled

end module bio_params