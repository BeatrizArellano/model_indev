! src/physics/physics_main.F90
module physics_main
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use load_forcing,     only: ForcingState, ForcingYearData, load_year_data, print_forcing_summary
  use tidal_parameters_readers, only: TidalParams, Constituent, read_tidal_parameters
  
  implicit none
  private


  ! Public API
  public :: init_physics, end_physics             ! set up grid, parameters, state

  !======================
  ! Internal state
  !======================


  type :: PhysState
     real(rk), allocatable :: T(:)      ! temperature [°C]
     !real(rk), allocatable :: S(:)      ! salinity [psu]
     real(rk), allocatable :: Kz(:)     ! vertical diffusivity [m2/s]
     ! add more prognostics/diagnostics as you grow
  end type

  ! Module-scope singletons (encapsulated)
  type(PhysState)    :: PState
  type(TidalParams)  :: Tides
  type(Constituent)  :: m2
  type(Constituent)  :: s2
  type(Constituent)  :: k1
  type(Constituent)  :: o1
  type(Constituent)  :: n2
  logical            :: is_phys_initialized = .false.

contains

  !======================
  ! Initialization
  !======================
  subroutine init_physics(cfg_params,location)
    ! Receives from main: user config and location
    type(ConfigParams), intent(in)    :: cfg_params
    type(LocationInfo), intent(in)    :: location    
    
    type(ForcingState)    :: FS
    type(ForcingYearData) :: surf

    character(len=:), allocatable :: filename, file_var    
    logical :: ok
    character(len=256) :: errmsg
    integer  :: i, y, k, nt, iu
    character(len=128) :: fname

    if (is_phys_initialized) return    
    ! -Reading  tidal parameters ----
    call read_tidal_parameters(cfg_params, location%lat,location%lon, Tides)   ! Reads tidal parameters
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
    
  


    ! ---- 3) Allocate & set initial state ----
    !call allocate_state(G, X)
    !call set_initial_conditions(cfg, G, X)

    

    ! ---- 5) Initialize turbulence/mixing scheme ----
    !call turbulence_init(cfg, G, P, X)

    is_phys_initialized = .true.
  end subroutine init_physics 


  !======================
  ! Finalize
  !======================
  subroutine end_physics()
    if (.not. is_phys_initialized) return
    !call deallocate_state(X)
    !call deallocate_grid(G)
    Tides = TidalParams()  ! lets allocatables inside tidy up via assignment
    is_phys_initialized = .false.
  end subroutine end_physics


  !subroutine dev_test()
  !    !--------------- Tests -----------------------
  !  
  !    do y = start_datetime%year, end_datetime%year
  !      k = y - FS%sim_y_start + 1
!
  !      call load_year_data(FS, k, surf, ok, errmsg)
  !      if (.not. ok) stop trim(errmsg)
!
  !      ! assume all vars same length and file-backed for this quick test
  !      nt = size(surf%air_temp%data)
!
  !      write(fname,'(A,I0,A)') 'forcing_', y, '.dat'   ! e.g., forcing_2001.dat
  !      open(newunit=iu, file=trim(fname), status='replace', action='write')
!
  !      write(iu,'(A)') 'surf_air_temp sl_pressure relative_humidity shortwave_radiation longwave_radiation wind_speed wind_direction co2_air'
  !      do i = 1, nt
  !        write(iu,'(8(1X,ES16.8))') surf%air_temp%data(i), surf%slp%data(i),       &
  !                                    surf%rel_hum%data(i),  surf%short_rad%data(i), &
  !                                    surf%long_rad%data(i), surf%wind_spd%data(i),  &
  !                                    surf%wind_dir%data(i), surf%co2_air%const_value
  !      end do
  !      close(iu)
  !    end do
!
  !    call load_year_data(FS, 1, surf, ok, errmsg)
  !    if (.not. ok) stop trim(errmsg)
!
  !    !call surf%air_temp%init_cursor()
  !    time = (30*86400)+10
  !    write(*,'(A,1X,ES12.5)') 'Tair_:',   surf%air_temp%value_at_step(time)
  !    time = (350*86400) + 86300
  !    write(*,'(A,1X,ES12.5)') 'Tair_180:', surf%air_temp%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'psl:', surf%slp%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'wind_speed:', surf%wind_spd%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'wind_dir:', surf%wind_dir%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'shortwave:', surf%short_rad%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'longwave:', surf%long_rad%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'rel_hum:', surf%rel_hum%value_at_step(time)
  !    write(*,'(A,1X,ES12.5)') 'co2:', surf%co2_air%value_at_step(time)
!
  !    
!
!
  !    filename       = cfg_params%get_param_str('forcing.filename',default='')
  !    ! Per-var filename: use it only if set (not null/missing/empty)
  !    if (.not. cfg_params%is_disabled('forcing.surf_air_temp.filename')) then
  !      file_var = cfg_params%get_param_str('forcing.surf_air_temp.filename', empty_ok=.false., trim_value=.true.)
  !    else
  !      file_var = filename
  !    end if
!
  !    write(*,*)  'filename=', adjustl(trim(file_var))
!
  !    if (.not. cfg_params%is_disabled('forcing.surf_air_temp.constant')) then     
  !      test_param = cfg_params%get_param_num('forcing.surf_air_temp.constant', finite=.true.)
  !      write(*,*)  'constant=', test_param
  !    else
  !      write(*,*)  'constant is Null and not active' 
  !    end if
  !end subroutine dev_test()

end module physics_main
