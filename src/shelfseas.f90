! src/physics/shelfseas.F90
module shelfseas
  use precision_types,  only: rk, lk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_types,       only: DateTime, CFCalendar
  use validation_utils, only: validate_input_dates, validate_location_input, print_header
  use sim_clocks,       only: init_clock
  use grids,            only: VerticalGrid, build_water_grid, write_vertical_grid
  use forcing_manager,  only: ForcingManager, ForcingSnapshot
  use physics_types,    only: PhysicsState, PhysicsEnv
  use physics_main,     only: init_physics, solve_physics, end_physics
  use output_manager,   only: OutputManager
  use bio_params,       only: is_bio_enabled
  use bio_types,        only: BioEnv
  use bio_main,         only: init_bio_fabm, integrate_bio_fabm, end_bio_fabm

  !Dev
  use time_utils, only:   datetime_to_str

  implicit none
  private

  !======================
  ! Public API
  !======================
  public :: init_shelfseas, run_shelfseas, end_shelfseas             ! set up grid, parameters, state 

  
  character(len=20), public :: config_file = 'main.yaml'  
  logical  :: bio_enabled         = .false.
  logical  :: sediments_enabled   = .false.
  logical  :: is_main_initialized = .false.
  logical  :: stop_on_error = .true.                 ! Hard stop on error  -> useful in the future when running multiple columns
  logical  :: ok = .false.
  character(len=256) :: msg

  integer(lk)   ::  dt, n_steps                      ! time-step (seconds) and total number of steps
  integer(lk)   ::  sim_length_sec, last_dt_length   ! Length of simulation in seconds and duration of last time-step

  type(ConfigParams)    :: cfg_params
  type(LocationInfo)    :: location
  type(DateTime)        :: start_datetime, end_datetime
  type(CFCalendar)      :: calendar
  type(VerticalGrid)    :: wgrid
  type(ForcingManager)  :: ForcMan
  type(ForcingSnapshot) :: ForcSnp
  type(PhysicsEnv)      :: PE
  type(BioEnv), target  :: BE
  type(OutputManager)   :: OM
  character(:), allocatable :: time_units, calname
  

contains

  !======================
  ! Initialization
  !======================

    subroutine init_shelfseas()
        integer, parameter :: cal_unknown = 0
        if (is_main_initialized) return
        ! Load user configuration in yaml file
        call cfg_params%init()    
        call cfg_params%load_yaml_content(config_file)
        ! Read and validate location parameters
        call validate_location_input(cfg_params, location)
        ! Read and validate simulation dates and calendar
        call validate_input_dates(cfg_params, start_datetime, end_datetime, calendar)        
        ! Print header for simulation
        call print_header(location,start_datetime,end_datetime) 

        call is_bio_enabled(cfg_params, bio_enabled, sediments_enabled)
        write(*,*) "Bio=", bio_enabled, ' Sed=', sediments_enabled

        ! Builds vertical grid
        call build_water_grid(cfg_params, location%depth, wgrid)
        call write_vertical_grid(wgrid, 'Vertical_grid.dat')
        ! Verify and initialise forcing data
        call ForcMan%set_hard_fail(stop_on_error)  

        call ForcMan%init(cfg_params, calendar, location, start_datetime, end_datetime, ok, msg)
        if (.not. ok) stop 'init_forcing failed: '//trim(msg)           
        ! Sets calendar from forcing data
        if (calendar%kind == cal_unknown) then
             calendar%kind = ForcMan%get_sim_calendar()
        end if
        ! Reference date is the start datetime for the simulation
        dt = cfg_params%get_param_int('time.dt', default=300, positive=.true.)        
        call init_clock(calendar,start_datetime, end_datetime, dt, sim_length_sec, n_steps, last_dt_length)
        write(*,*) 'n_steps=', n_steps, 'last=', last_dt_length
        call ForcMan%prepare(dt, preload_pad_sec=3600_lk, ok=ok, errmsg=msg)
        if (.not. ok) stop 'prepare_forcing failed: '//trim(msg)
        ! Initialise physics
        call init_physics(cfg_params, location, wgrid, PE)

        if (bio_enabled) then
            call init_bio_fabm(cfg_params, location, wgrid, dt, PE%PS, ForcSnp, BE)
        end if
        !!!! In dev
        time_units = 'seconds since ' // trim(datetime_to_str(start_datetime))
        calname    = trim(calendar%name())
        call OM%init(cfg_params, PE%grid, dt_s=dt, time_units=time_units, calendar_name=calname, &
                     title='ShelfSeas output')

        is_main_initialized = .true.     
        
    end subroutine init_shelfseas

    ! Main subroutine to run ShelfSeas. 
    ! Main time loop is here. 
    subroutine run_shelfseas()
        implicit none

        integer(lk)        :: istep, elapsed_time
        integer            :: ierr
        character(len=512) :: errmsg
        logical            :: is_first_step = .true.  

        if (.not. is_main_initialized) error stop 'run_shelfseas: shelfseas not initialised.'

        ! Main simulation loop
        do istep = 1_lk, n_steps
            if (istep > 1_lk) is_first_step = .false.
            elapsed_time = (istep - 1_lk) * dt
            call ForcMan%tick(elapsed_time)                                   ! Time-manager to load yearly forcing data on time
            call ForcMan%sample(elapsed_time, ForcSnp)                        ! get forcing snapshot for the current model time
            call solve_physics(PE, ForcSnp, dt, elapsed_time, is_first_step, ierr, errmsg)
            if (ierr /= 0) then
                write(*,*) 'solve_physics failed: ', trim(errmsg)
                return
            end if
            if (bio_enabled) then
                call integrate_bio_fabm(BE, PE%PS, ForcSnp, dt, istep)
            end if
            call OM%step(PE%PS%temp, PE%PS%Kz) 
        end do
    end subroutine run_shelfseas

    subroutine end_shelfseas()
        if (.not. is_main_initialized) return
        call OM%close(sync_now=.true.)
        call end_physics(PE)
        if (bio_enabled) then
            call end_bio_fabm(BE)
        end if
        call cfg_params%clear()
    end subroutine end_shelfseas 

end module shelfseas
