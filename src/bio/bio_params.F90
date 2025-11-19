module bio_params
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams

  implicit none
  private

  public :: is_bio_enabled

type, public :: BioParams
    logical :: sediments_enabled = .false.
end type BioParams


contains
    
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