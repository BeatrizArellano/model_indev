! src/physics/physics_driver.F90
module physics_driver
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_utils,       only: DateTime
  use calendar_types,   only: CFCalendar
  use tidal_parameters_readers, only: TidalParams, Constituent, read_tidal_parameters
  
  implicit none
  private

  !======================
  ! Public API
  !======================
  public :: physics_init             ! set up grid, parameters, state
  public :: physics_end              ! release resources
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

  type :: PhysState
     real(rk), allocatable :: T(:)      ! temperature [°C]
     !real(rk), allocatable :: S(:)      ! salinity [psu]
     real(rk), allocatable :: Kz(:)     ! vertical diffusivity [m2/s]
     ! add more prognostics/diagnostics as you grow
  end type

  ! Module-scope singletons (encapsulated)
  type(Grid1D)       :: Grid
  type(PhysState)    :: PState
  type(TidalParams)  :: Tides
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
  subroutine physics_init(location, start_datetime, end_datetime, calendar)
    ! Receives from main: user config, simulation dates and calendar   
    type(LocationInfo), intent(in) :: location
    type(DateTime),     intent(in) :: start_datetime, end_datetime
    type(CFCalendar),   intent(in) :: calendar    

    type(ConfigParams) :: physics_cfg


    character(len=:), allocatable :: filename, file_var
    real(rk) :: test_param
    integer  :: i

    call physics_cfg%init()    
    call physics_cfg%load_yaml_content('physics.yaml')

    filename       = physics_cfg%get_param_str('forcing.filename',default='')
    ! Per-var filename: use it only if set (not null/missing/empty)
    if (physics_cfg%is_set('forcing.surf_air_temp.filename')) then
      file_var = physics_cfg%get_param_str('forcing.surf_air_temp.filename', empty_ok=.false., trim_value=.true.)
    else
      file_var = filename
    end if

    write(*,*)  'filename=', adjustl(trim(file_var))


    if (.not. physics_cfg%is_null('forcing.surf_air_temp.constant')) then     
      test_param = physics_cfg%get_param_num('forcing.surf_air_temp.constant', finite=.true.)
      write(*,*)  'constant=', test_param
    else
      write(*,*)  'constant is Null and not active'
    end if



    if (is_initialized) return

    
    ! ---- 2) Build grid ----
    !call build_vertical_grid(nz, P%h, G)

    ! ---- 3) Allocate & set initial state ----
    !call allocate_state(G, X)
    !call set_initial_conditions(cfg, G, X)

    ! ---- 4) Read tidal parameters once ----
    call read_tidal_parameters(physics_cfg, location%lat,location%lon, Tides)   ! Reads tidal parameters
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
