module bio_types
    use bio_params,        only: BioParams, SedParams
    use event_manager,     only: EventManager   
    use fabm,              only: type_fabm_model, type_fabm_interior_variable_id, &
                                 type_fabm_horizontal_variable_id, type_fabm_scalar_variable_id
    use grids,             only: VerticalGrid
    use precision_types,   only: rk, lk
    use tridiagonal,       only: TridiagCoeff    
    use variable_registry, only: VarMetadata
    

  implicit none

  integer, parameter, public :: DIFF_NONE=0, DIFF_O2CO2_AB=1, DIFF_ION_LINEAR=2, &
                                DIFF_ARRHENIUS=3, DIFF_WILKE_CHANG=4, DIFF_STOKES_EINSTEIN=5

  public :: BioState, BioEnv, SedimentEnv

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
    real(rk), allocatable :: temp(:), sal(:), rho(:)         ! temperature, salinity and density (From bottom to surface)
    real(rk), allocatable :: pres(:)                         ! pressure [dbar]     (optional)
    real(rk), allocatable :: swr(:),  par(:)                 ! PAR profile [W/m2]  (optional)      
    real(rk), allocatable :: atten_coeff(:)                  ! [m-1] attenuation coefficient for PAR, size N
    real(rk), allocatable :: vert_diff(:)                    ! vertical diffusivity
    ! Surface variables
    real(rk) :: short_rad  = 0._rk
    real(rk) :: par_sfc    = 0._rk                          ! surface PAR [W/m2] (optional)
    real(rk) :: wind_spd   = 0._rk
    real(rk) :: co2_air    = 0._rk
    real(rk) :: slp        = 0._rk
    real(rk) :: cloud      = 0._rk                          ! cloud fraction [0-1], optional
    real(rk) :: stressb    = 0._rk                          ! bottom stress magnitude [Pa], optional
    real(rk) :: ice_af     = 0._rk                          ! Ice area fraction (Always 0 for now)
    ! Bottom variables for the sediments
    real(rk) :: u_taub     = 0._rk                          ! Friction velocity
    real(rk) :: z0b        = 0._rk                          ! Bottom roughness length
    real(rk) :: Nz_btm     = 0._rk                          ! momentum viscosity just above the bottom
    real(rk) :: Kz_btm     = 0._rk                          ! Eddy diffusivity just above the bottom
    ! Number of days since the start of the year
    real(rk) :: doy        = 0._rk
  end type BioState

  type, public :: TracerProperties
      integer  :: fabm_index      = -1          ! index in interior_state_variables
      logical  :: is_solute       = .false.
      logical  :: is_particulate  = .false.
      logical  :: disable_transport = .false.
      ! Adsorption 
      real(rk) :: adsorption      = 0._rk       ! Coefficient for adsorption on sediment grains

      !--- Diffusivity properties retrieved from FABM
      integer  :: diff_method = DIFF_NONE       ! Selects how molecular diffusivity is computed for solute tracers.
      ! O2 / CO2 empirical formulation (Boudreau 1997, Eqs. 4.58–4.59)
      real(rk) :: A = 0._rk                     ! 
      real(rk) :: B = 0._rk                     ! Slope of the linear relationship
      ! Linear regresions of the Infinite-Dilution Diffusion Coefficients Do for ions vs Temperature       
      ! Coefficients m0, m1 are taken directly from Boudreau (1997) Tables 4.7 and 4.8. 
      !-- Ions: D = m0 + m1*t
      real(rk) :: m0 = 0._rk
      real(rk) :: m1 = 0._rk      
      !-- Arrhenius: D = A0 * exp(-Ea/(R*T)) (Eq. 4.60 in Boudreau,1997)
      real(rk) :: A0 = 0._rk                   ! units 10-5 cm2 s-1   (Following the units in Boudreau Table 4.4)
      real(rk) :: Ea = 0._rk                   ! kJ mol-1
      ! Wilke–Chang (Wilke & Chang 1955; Hayduk & Laudie 1974 modification)
      ! 4.72E-09 * TK / (mu * Vb^0.6)  (Eq. 4.57 in Boudreau)
      ! Numerical prefactor and cm2->m2 conversion applied in compute routine.
      real(rk) :: Vb = 0._rk                   ! molar volume at boiling point [cm3/mol] (typical)
      ! Stokes–Einstein scaling from a reference value (Sref=0, Pref=1 atm):
      ! D(T,S,P) = Dref * (TK/TrefK) * (mu(Tref,0,Patm) / mu(T,S,P))
      real(rk) :: Dref = 0._rk                 ! Reference diffusivity [m^2 s^-1] measured in pure water (S=0) at atmospheric pressure
      real(rk) :: Tref = 0._rk                 ! Reference temperature [°C]    
      real(rk) :: Sref = 0._rk                 ! Reference salinity [PSU]     
   end type TracerProperties


  !======================
  ! Sediment Environment
  !======================
  type, public :: SedimentEnv
      logical :: is_init = .false.

      ! A pointer to the sediment grid owned by BioEnv
      type(VerticalGrid), pointer :: grid => null()

      ! Parameters
      ! IMPORTANT:  Internal sediment computations must use params_SI only.
      ! params_user stores (only for reference) the parameters in units commonly reported in the literature.
      type(SedParams) :: params_user       ! Sediment configuation parameters in commonly reported units (cm, yr)
      type(SedParams) :: params_SI         ! Sediment configuation parameters in standard units (m and s)      

      integer :: nz = 0
      ! ---- Properties at layers' centres (1:nz)
      real(rk), allocatable :: poro(:)              ! [-] porosity at centres
      real(rk), allocatable :: porewat_thickness(:) ! [m] Storage capacity of porewater per unit horizontal-area
      real(rk), allocatable :: solid_thickness(:)   ! [m] Storage capacity of solids per unit horizontal area
      real(rk), allocatable :: bioirr(:)            ! [s-1] 

      ! ---- Properties at layer interfaces (0:nz)
      real(rk), allocatable :: poro_w(:)     ! [-] porosity at interfaces
      real(rk), allocatable :: theta2(:)      ! [-] Diffusion tortuosity factor (Boudreau 1997), used as D_eff = D0 / theta2
      real(rk), allocatable :: bioturb(:)    ! [m2/s] particulate diffusivity    
      real(rk), allocatable :: bioirr_w(:)   ! [s-1] 
      ! --- Burial velocities
      real(rk), allocatable :: vel_solids(:)     ! [m/s] Burial velocity for particulate matter
      real(rk), allocatable :: vel_solutes(:)    ! [m/s] Burial velocity for porewater

      ! ---- Working arrays
      real(rk), allocatable :: bulk_conc(:,:)      ! (nsed, n_interior) Bulk-sediment concentrations for mass-conservation
      real(rk), allocatable :: diff_sed(:)         ! Array to store effective diffusivities scaled by tortuosity per tracer
      real(rk), allocatable :: diff_sed0(:)        ! Array to store free-solution diffusivities (not scaled)
      real(rk), allocatable :: diff_sed_max(:)     ! Array to store maximum effective diffusivities scaled by tortuosity for all tracers
      real(rk), allocatable :: Db_eff_solids(:)    ! Array to store effective bioturbation diffusivity (scaled by 1-phi)
      real(rk), allocatable :: swi_flux(:)         ! Flux of solutes at the sediment-water interface (tracer specific)
      ! Working space
    type(TridiagCoeff)      :: sed_trid            ! workspace for solving scalar diffusion
  end type SedimentEnv  

  ! An envelope for the Biogeochemistry Environment in one column
  type, public :: BioEnv
    class (type_fabm_model), pointer :: model => null()        ! FABM model instance
    !----------- State and parameters
    type(BioState)         :: BS                               ! State of tracers
    type(BioParams)        :: params                           ! Biogeochemical configuration parameters
    ! --- Sediment environment
    type(SedimentEnv)      :: SED
    ! ---------- Grids
    type(VerticalGrid)     :: grid                             ! full column grid 
    type(VerticalGrid)     :: wat_grid                         ! water-column grid 
    type(VerticalGrid)     :: sed_grid                         ! sediments grid 
    integer :: nsed = 0                                        ! Number of layers in sediments
    integer :: nwat = 0                                        ! Number of layers in water
    integer :: k_sed_btm = 0, k_sed_sfc = 0                    ! Bottom and surface Indices for sediments
    integer :: k_wat_btm = 0, k_wat_sfc = 0                    ! Bottom and surface Indices for water
    !---------- Event Manager
    type(EventManager) :: Events

    !------------ FABM environment variable ids
    type (type_fabm_interior_variable_id)   :: id_temp, id_salt, id_rho, id_swr, id_par, id_pres, id_atten
    type (type_fabm_horizontal_variable_id) :: id_windspd, id_par_sfc, id_slp, id_cloud, id_stressb, id_swr_sfc, id_co2, id_ice_af
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
    logical :: need_ice_af  = .false.                          ! Ice area fraction
    logical :: is_init = .false.                               ! has biogeochemistry been initialised?
    ! Pointer to attenuation coefficient
    real(rk), pointer :: atten_ptr(:) => null()
    ! Working space
    type(TridiagCoeff)     :: wat_trid                             ! workspace for solving scalar diffusion
    ! Interior diagnostics 
    integer :: n_diag_int = 0
    integer,   allocatable :: diag_int_index(:)                ! FABM indices of interior diagnostics
    real(rk),  allocatable :: diag_int(:,:)                    ! Array for interior diagnostics(nz, n_diag_int)
    type(VarMetadata), allocatable :: diag_int_vars(:)         ! Metadata for interior diagnostic variables
    type(TracerProperties), allocatable :: tracer_info(:)
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
    ! Other environment arrays
    type(VarMetadata), allocatable :: env_int_vars(:)          ! 1D environment variables like par/swr/pres, etc.
    ! Working arrays
    real(rk), allocatable :: velocity(:,:)                     ! Vertical velocity due to residual movement (nz, n_interior)
    real(rk), allocatable :: vel_faces(:,:)                    ! Vertical velocity due to residual movement at interfaces (0:nz, n_interior)
    real(rk), allocatable :: tendency_int(:,:)                 ! Tendencies (source terms) for interior tracers from FABM
    real(rk), allocatable :: tendency_sf(:), tendency_bt(:)    ! source terms at the surface and bottom. They store the tendencies (or derivatives dC/dt) for the tracers
    real(rk), allocatable :: tendency_int_evt(:,:)             ! Tendencies for interior tracers resulting from external events (not from FABM)
    real(rk), allocatable :: tendency_int_total(:,:)           ! Total tendencies for interior tracers from FABM and from external events 
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