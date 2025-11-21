module bio_types
    use precision_types, only: rk, lk
    use grids,           only: VerticalGrid
    use fabm,            only: type_fabm_model, type_fabm_interior_variable_id, &
                                   type_fabm_horizontal_variable_id, type_fabm_scalar_variable_id
    use tridiagonal,     only: TridiagCoeff
    use bio_params,      only: BioParams
    

  implicit none

  public :: BioState, BioEnv

  !======================
  ! Internal state
  !======================

  type, public :: BioState
    integer :: n_interior = 0                                ! Number of interior biogeochemical variables
    integer :: n_surface  = 0                                ! Number of surface only variables
    integer :: n_bottom   = 0                                ! Number of bottom only variables
    integer :: n_total    = 0                                ! Total number of variables

    character(len=:), allocatable :: intvar_names(:)         ! Names for interior biogeochemical variables
    character(len=:), allocatable :: sfcvar_names(:)         ! Names for surface biogeochemical variables
    character(len=:), allocatable :: btmvar_names(:)         ! Names for bottom biogeochemical variables
    ! State arrays
    real(rk),         allocatable :: interior_state(:,:)     ! State of interior variables (nz,n_interior)
    real(rk),         allocatable :: bottom_state(:)         ! State of bottom variables
    real(rk),         allocatable :: surface_state(:)        ! State of surface variables
    ! Relevant Physical variables
    real(rk), allocatable :: temp(:), sal(:), rho(:)         ! temperature, salinity and density (From surface to bottom)
    real(rk), allocatable :: pres(:)                         ! pressure [dbar]       (optional)
    real(rk), allocatable :: swr(:),  par(:)                 ! PAR profile [W/m2]  (optional)
    real(rk), allocatable :: vert_diff(:)                    ! vertical diffusivity
    ! Surface variables
    real(rk) :: short_rad  = 0._rk
    real(rk) :: par_sfc    = 0._rk                          ! surface PAR [W/m2] (optional)
    real(rk) :: wind_spd   = 0._rk
    real(rk) :: co2_air    = 0._rk
    real(rk) :: slp        = 0._rk
    real(rk) :: cloud      = 0._rk                          ! cloud fraction [0-1], optional
    real(rk) :: stressb    = 0._rk                          ! bottom stress magnitude [Pa], optional
    ! Number of days since the start of the year
    real(rk) :: doy        = 0._rk
  end type BioState

  ! An envelope for the Biogeochemistry Environment in one column
  type, public :: BioEnv
    class (type_fabm_model), pointer :: model => null()  ! FABM model instance
    type(VerticalGrid)     :: grid                             ! full column grid 
    type(VerticalGrid)     :: wat_grid                             ! water column grid 
    type(VerticalGrid)     :: sed_grid                             ! sediments grid 
    type(BioState)         :: BS                               ! State of tracers
    type(BioParams)        :: params                           ! per-column biogeochemical params (or shared elsewhere)
    ! FABM environment variable ids
    type (type_fabm_interior_variable_id)   :: id_temp, id_salt, id_rho, id_swr, id_par, id_pres
    type (type_fabm_horizontal_variable_id) :: id_windspd, id_par_sfc, id_slp, id_cloud, id_stressb, id_swr_sfc, id_co2
    type (type_fabm_scalar_variable_id)     :: id_yearday
    ! Which env vars are actually needed?
    logical :: need_temp    = .false.
    logical :: need_salt    = .false.
    logical :: need_rho     = .false.
    logical :: need_pres    = .false.
    logical :: need_par     = .false.
    logical :: need_swr     = .false.
    logical :: need_windspd = .false.
    logical :: need_slp     = .false.
    logical :: need_par_sfc = .false.
    logical :: need_swr_sfc = .false.
    logical :: need_cloud   = .false.
    logical :: need_stressb = .false.
    logical :: need_co2     = .false.
    logical :: is_init = .false.                               ! has biogeochemistry been initialised?
    ! Working space
    type(TridiagCoeff)     :: trid                             ! workspace for implicit solves in scalar diffusion
    ! Working arrays
    real(rk), allocatable :: velocity(:,:)                     ! Vertical velocity due to residual movement (nz, n_interior)
    real(rk), allocatable :: tendency_int(:,:)                 ! Tendencies (source terms) for interior tracers
    real(rk), allocatable :: tendency_sf(:), tendency_bt(:)    ! source terms at the surface and bottom. They store the tendencies (or derivatives dC/dt) for the tracers
    real(rk), allocatable :: flux_sf(:), flux_bt(:)            ! Fluxes of interior variables at the surface and bottom
    ! Counters for how many times the state was repaired
    integer :: nrepair_int = 0
    integer :: nrepair_sfc = 0
    integer :: nrepair_btm = 0
    logical :: valid_int = .true.
    logical :: valid_sfc = .true.
    logical :: valid_btm = .true.
  end type BioEnv


end module bio_types