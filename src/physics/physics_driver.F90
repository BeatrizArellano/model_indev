! src/physics/physics_driver.F90
module physics_driver
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_utils,       only: DateTime
  use calendar_types,   only: CFCalendar
  use load_forcing,     only: ForcingState, ForcingYearData, init_forcing, load_year_data, print_forcing_summary
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
    type(ForcingState) :: FS
    type(ForcingYearData):: surf



    character(len=:), allocatable :: filename, file_var
    real(rk) :: test_param, time
    logical :: ok
    character(len=256) :: errmsg
    integer  :: i, y, k, nt, iu
    character(len=128) :: fname

    if (is_initialized) return

    call physics_cfg%init()    
    call physics_cfg%load_yaml_content('physics.yaml')

    ! ------------------ Forcing data -------------------------------------
    write(*,'(A)') 'Scanning forcing data...'
    call init_forcing(physics_cfg, calendar, location, start_datetime, end_datetime, FS, ok, errmsg)
    if (.not. ok) then
      write(*,*) trim(errmsg)  
      stop                    
    end if
    call print_forcing_summary(FS)
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
    
    !--------------- Tests -----------------------
    
    do y = start_datetime%year, end_datetime%year
      k = y - FS%sim_y_start + 1

      call load_year_data(FS, k, surf, ok, errmsg)
      if (.not. ok) stop trim(errmsg)

      ! assume all vars same length and file-backed for this quick test
      nt = size(surf%air_temp%data)

      write(fname,'(A,I0,A)') 'forcing_', y, '.dat'   ! e.g., forcing_2001.dat
      open(newunit=iu, file=trim(fname), status='replace', action='write')

      write(iu,'(A)') 'surf_air_temp sl_pressure relative_humidity shortwave_radiation longwave_radiation wind_speed wind_direction co2_air'
      do i = 1, nt
        write(iu,'(8(1X,ES16.8))') surf%air_temp%data(i), surf%slp%data(i),       &
                                    surf%rel_hum%data(i),  surf%short_rad%data(i), &
                                    surf%long_rad%data(i), surf%wind_spd%data(i),  &
                                    surf%wind_dir%data(i), surf%co2_air%const_value
      end do
      close(iu)
    end do

    call load_year_data(FS, 1, surf, ok, errmsg)
    if (.not. ok) stop trim(errmsg)

    !call surf%air_temp%init_cursor()
    time = (30*86400)+10
    write(*,'(A,1X,ES12.5)') 'Tair_:',   surf%air_temp%value_at_step(time)
    time = (359*86400) + 86300
    write(*,'(A,1X,ES12.5)') 'Tair_180:', surf%air_temp%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'psl:', surf%slp%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'wind_speed:', surf%wind_spd%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'wind_dir:', surf%wind_dir%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'shortwave:', surf%short_rad%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'longwave:', surf%long_rad%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'rel_hum:', surf%rel_hum%value_at_step(time)
    write(*,'(A,1X,ES12.5)') 'co2:', surf%co2_air%value_at_step(time)

    


    filename       = physics_cfg%get_param_str('forcing.filename',default='')
    ! Per-var filename: use it only if set (not null/missing/empty)
    if (.not. physics_cfg%is_disabled('forcing.surf_air_temp.filename')) then
      file_var = physics_cfg%get_param_str('forcing.surf_air_temp.filename', empty_ok=.false., trim_value=.true.)
    else
      file_var = filename
    end if

    write(*,*)  'filename=', adjustl(trim(file_var))

    if (.not. physics_cfg%is_disabled('forcing.surf_air_temp.constant')) then     
      test_param = physics_cfg%get_param_num('forcing.surf_air_temp.constant', finite=.true.)
      write(*,*)  'constant=', test_param
    else
      write(*,*)  'constant is Null and not active'
    end if



    

    
    ! ---- 2) Build grid ----
    !call build_vertical_grid(nz, P%h, G)

    ! ---- 3) Allocate & set initial state ----
    !call allocate_state(G, X)
    !call set_initial_conditions(cfg, G, X)

    

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
