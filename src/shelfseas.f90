! src/physics/shelfseas.F90
module shelfseas
  use precision_types,  only: rk, lk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_types,       only: DateTime, CFCalendar
  use validation_utils, only: validate_input_dates, validate_location_input, print_header
  use sim_clocks,       only: init_clock
  use grid_builders,    only: VerticalGrid, build_water_grid, write_vertical_grid
  use forcing_manager, only: ForcingManager, ForcingSnapshot
  use physics_main,     only: init_physics, end_physics
  
  implicit none
  private

  !======================
  ! Public API
  !======================
  public :: init_shelfseas, run_shelfseas, end_shelfseas             ! set up grid, parameters, state 

  
  character(len=20), public :: config_file = 'main.yaml'  
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
        call init_physics(cfg_params,location)
        is_main_initialized = .true.
    end subroutine init_shelfseas

    ! Main subroutine to run ShelfSeas. 
    ! Main time loop is here. 
    subroutine run_shelfseas()
        implicit none
        integer(lk) :: i, t_step_sec, dt_win

        !----------- For testing--------------------------------------------
        integer, parameter :: SECS_PER_DAY = 86400
        integer, parameter :: n_tests = 6
        integer, parameter :: test_days(n_tests) = [1, 7, 157, 360, 361, 500]
        logical            :: printed(n_tests)
        integer            :: k, day_num
        printed = .false.
        !---------------------------------------------------------------------


        if (.not. is_main_initialized) error stop 'run_shelfseas: shelfseas not initialised.'
        do i = 1_lk, n_steps
            ! Model time at the start of this step (seconds since start_datetime)
            t_step_sec = (i - 1_lk) * dt
            call ForcMan%tick(t_step_sec)                    ! Time-manager to load yearly forcing data on time
            call ForcMan%sample(t_step_sec, ForcSnp)         ! get forcing snapshot for the current model time

            !-------------- TEST --------------------------------------------------------
            ! Determine (1-based) day number for current step
            day_num = int(t_step_sec / SECS_PER_DAY, kind=lk) + 1

            ! If this is one of the target days and we haven't printed yet, dump the snapshot
            do k = 1, n_tests
            if (.not. printed(k) .and. day_num == test_days(k)) then
                write(*,'(A,I0)') '--- Forcing snapshot at day ', test_days(k)
                write(*,'(A,F12.5)') '  air_temp       = ', ForcSnp%air_temp
                write(*,'(A,F12.5)') '  slp            = ', ForcSnp%slp
                write(*,'(A,F12.5)') '  rel_hum        = ', ForcSnp%rel_hum
                write(*,'(A,F12.5)') '  short_rad      = ', ForcSnp%short_rad
                write(*,'(A,F12.5)') '  long_rad       = ', ForcSnp%long_rad
                write(*,'(A,F12.5)') '  wind_spd       = ', ForcSnp%wind_spd
                write(*,'(A,F12.5)') '  wind_dir       = ', ForcSnp%wind_dir
                write(*,'(A,F12.5)') '  co2_air        = ', ForcSnp%co2_air
                printed(k) = .true.
            end if
            end do
            !--------------------------------------------------------------------------

            ! Check Forcing: ensure active data covers [t_step_sec, t_step_sec + dt_win]            

            ! Physics step (handles its own subcycling/implicit solves as needed)
            ! call physics_step(t_step_sec, dt_win, wgrid, ...)

            ! Transport step (only if enabled; also subcycles internally)
            ! call transport_step(t_step_sec, dt_win, wgrid, ...)

            ! Output/diagnostics if due
            ! call output_maybe_write(t_step_sec, dt_win, ...)
        end do


    end subroutine run_shelfseas

    subroutine end_shelfseas()
        if (.not. is_main_initialized) return
        call end_physics()
        call cfg_params%clear()
    end subroutine end_shelfseas 

end module shelfseas
