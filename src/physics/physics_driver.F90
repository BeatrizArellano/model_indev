! src/physics/physics_driver.F90
module physics_driver
  use precision_types, only: rk
  use read_config_yaml, only: ConfigParams
  use tidal_parameters_readers, only: TidalParams, Constituent, read_tidal_parameters
  implicit none
  private

  !======================
  ! Public API
  !======================
  public :: physics_init            ! set up grid, parameters, state
  public :: physics_end        ! release resources
  !public :: physics_get_grid        ! optional accessor(s)

  !======================
  ! Internal state
  !======================
  type :: Grid1D
     integer :: nz = 0
     real(rk), allocatable :: z(:)     ! cell centers [m]
     real(rk), allocatable :: dz(:)    ! layer thickness [m]
     real(rk), allocatable :: z_w(:)   ! interfaces (optional)
  end type

  type :: LocParams
     real(rk) :: lat
     real(rk) :: lon
     real(rk) :: depth      ! depth [m]
  end type

  type :: PhysState
     real(rk), allocatable :: T(:)      ! temperature [°C]
     !real(rk), allocatable :: S(:)      ! salinity [psu]
     real(rk), allocatable :: Kz(:)     ! vertical diffusivity [m2/s]
     ! add more prognostics/diagnostics as you grow
  end type

  ! Module-scope singletons (encapsulated)
  type(Grid1D)     :: Grid
  type(LocParams)  :: Loc
  type(PhysState)  :: PState
  type(TidalParams):: Tides
  type(Constituent)  :: m2
  type(Constituent)  :: s2
  type(Constituent)  :: k1
  type(Constituent)  :: o1
  type(Constituent)  :: n2
  logical          :: is_initialized = .false.

contains

  !======================
  ! Initialization
  !======================
  subroutine physics_init(user_cfg)    
    type(ConfigParams), intent(in) :: user_cfg
    real(rk) :: latitude, longitude, depth
    integer :: i


    if (is_initialized) return

    ! ---- 1) Read site/config ----
    latitude = user_cfg%get_param_num('main.location.latitude')
    longitude = user_cfg%get_param_num('main.location.longitude')    
    depth = user_cfg%get_param_num('main.location.depth') ! read depth as a parameter

    Loc%lat  = latitude
    Loc%lon  = longitude
    Loc%depth = depth

    write(*,*) 'latitude=', latitude, ' longitude=', longitude, ' depth=', depth

    ! ---- 2) Build grid ----
    !call build_vertical_grid(nz, P%h, G)

    ! ---- 3) Allocate & set initial state ----
    !call allocate_state(G, X)
    !call set_initial_conditions(cfg, G, X)

    ! ---- 4) Read tidal parameters once ----
    call read_tidal_parameters(user_cfg, Tides)   ! Reads tidal parameters
    call Tides%get('m2', m2)
    call Tides%get('s2', s2)
    call Tides%get('k1', k1)
    call Tides%get('o1', o1)
    call Tides%get('n2', n2)
    print *, '-----------------------------------------------'
    print *, 'Tidal constituents loaded: ', size(Tides%c)
    print *, '-----------------------------------------------'
    do i = 1, size(Tides%c)
      write(*,'(A3, 2X, "SEMA=",F8.3, 2X, "SEMI=",F8.3, 2X, &
              "INC=",F7.2," deg", 2X, "PHA=",F7.2," deg")') &
          trim(Tides%c(i)%name), Tides%c(i)%sema, Tides%c(i)%semi, &
          Tides%c(i)%inc_deg, Tides%c(i)%pha_deg
    end do
    print *, '-----------------------------------------------'



    ! ---- 5) Initialize turbulence/mixing scheme ----
    !call turbulence_init(cfg, G, P, X)

    is_initialized = .true.
  end subroutine physics_init

 
!
  !======================
  ! Finalize
  !======================
  subroutine physics_end()
    if (.not. is_initialized) return
    !call deallocate_state(X)
    !call deallocate_grid(G)
    Tides = TidalParams()  ! lets allocatables inside tidy up via assignment
    is_initialized = .false.
  end subroutine physics_end

end module physics_driver
