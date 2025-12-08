module bio_types
    use precision_types,   only: rk, lk
    use grids,             only: VerticalGrid
    use fabm,              only: type_fabm_model, type_fabm_interior_variable_id, &
                                   type_fabm_horizontal_variable_id, type_fabm_scalar_variable_id
    use tridiagonal,       only: TridiagCoeff
    use bio_params,        only: BioParams
    use variable_registry, only: VarMetadata
    

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
    type(VerticalGrid)     :: wat_grid                         ! water column grid 
    type(VerticalGrid)     :: sed_grid                         ! sediments grid 
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
    type(TridiagCoeff)     :: trid                             ! workspace for solving scalar diffusion
    ! Interior diagnostics 
    integer :: n_diag_int = 0
    integer,   allocatable :: diag_int_index(:)                ! FABM indices of interior diagnostics
    real(rk),  allocatable :: diag_int(:,:)                    ! Array for interior diagnostics(nz, n_diag_int)
    type(VarMetadata), allocatable :: diag_int_vars(:)         ! Metadata for interior diagnostic variables
    ! Horizontal diagnostics 
    integer :: n_diag_hz  = 0
    integer,   allocatable :: diag_hz_index(:)                 ! FABM indices of horizontal diagnostics
    real(rk),  allocatable :: diag_hz(:)                       ! Array for horizontal diagnostics (n_hdiag)
    type(VarMetadata), allocatable :: diag_hz_vars(:)          ! Metadata for horizontal diagnostic variables
    type(VarMetadata), allocatable :: int_vars(:)              ! Metadata for Biogeochemical interior variables
    type(VarMetadata), allocatable :: sfc_vars(:)              ! Metadata for surfcace variables
    type(VarMetadata), allocatable :: btm_vars(:)              ! Metadata for bottom variables
    ! conserved quantities diagnostics
    integer :: n_conserved = 0                                 ! number of FABM conserved quantities
    real(rk), allocatable :: conserved_interior(:,:)           ! (nz, n_conserved): interior values [per m3]
    real(rk), allocatable :: conserved_boundary(:)             ! (n_conserved): surface+bottom totals [per m2]
    real(rk), allocatable :: conserved_total(:)                ! column-integrated totals [per m2]
    type(VarMetadata), allocatable :: conserved_vars(:)        ! metadata for conserved totals
    ! Working arrays
    real(rk), allocatable :: velocity(:,:)                     ! Vertical velocity due to residual movement (nz, n_interior)
    real(rk), allocatable :: vel_faces(:,:)                    ! Vertical velocity due to residual movement at interfaces (0:nz, n_interior)
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