! src/physics/shelfseas.F90
module shelfseas
  use bio_main,          only: init_bio_fabm, integrate_bio_fabm, end_bio_fabm
  use bio_types,         only: BioEnv
  use event_manager,     only: EventManager   
  use geo_utils,         only: LocationInfo
  use grids,             only: VerticalGrid, build_grids
  use output_manager,    only: OutputManager
  use output_static,     only: StaticProfile, clear_static_profiles
  use physics_forcing,   only: PhysicsForcing, ForcingSnapshot
  use physics_main,      only: init_physics, solve_physics, end_physics
  use physics_types,     only: PhysicsState, PhysicsEnv
  use precision_types,   only: rk, lk
  use read_config_yaml,  only: ConfigParams
  use sim_clocks,        only: init_clock, print_progress, simtime_to_datetime 
  use state_loader,      only: StateData, load_state_file
  use time_types,        only: DateTime, CFCalendar, cal_unknown
  use time_utils,        only: datetime_to_str
  use validation_utils,  only: validate_input_dates, validate_location_input, print_header
  use variable_registry, only: VarMetadata 

  implicit none
  private

  !======================
  ! Public API
  !======================
  public :: init_shelfseas, run_shelfseas, end_shelfseas             ! set up grid, parameters, state 

  
  character(len=20), public :: config_file = 'main.yaml'  
  logical  :: is_bio_enabled      = .false.
  logical  :: is_sed_enabled      = .false.
  logical  :: is_main_initialized = .false.
  logical  :: stop_on_error = .true.                 ! Full stop on error  -> useful in the future when running multiple columns
  logical  :: load_yearly   = .false.                ! To load full series at the beginning or load every year
  logical  :: load_init_state = .false.              ! Whether Initial profiles need to be loaded to e.g. restart a simulation
  logical  :: ok = .false.
  character(len=256) :: msg

  integer(lk)   ::  dt, n_steps                      ! time-step (seconds) and total number of steps
  integer(lk)   ::  sim_length_sec, last_dt_length   ! Length of simulation in seconds and duration of last time-step
  integer :: nsed

  type(ConfigParams)    :: cfg_params
  type(LocationInfo)    :: location
  type(DateTime)        :: start_datetime, end_datetime, current_datetime
  type(CFCalendar)      :: calendar
  type(VerticalGrid)    :: wat_grid, sed_grid, full_grid
  type(PhysicsForcing)  :: PhysForc
  type(ForcingSnapshot) :: ForcSnp
  type(PhysicsEnv)      :: PE
  type(BioEnv), target  :: BE
  type(EventManager)    :: EVT 
  type(OutputManager)   :: OM
  type(StateData)       :: init_state
  type(VarMetadata),   allocatable :: all_vars(:)
  type(StaticProfile), allocatable :: static_profs(:)
  character(:),        allocatable :: time_units, calname
  character(len=:),    allocatable :: init_state_path     
  

contains

    !======================
    ! Initialisation
    !======================
    subroutine init_shelfseas()
        if (is_main_initialized) return
        ! Load user configuration in yaml file
        call cfg_params%init()    
        call cfg_params%load_yaml_content(config_file)
        ! Read and validate location parameters
        call validate_location_input(cfg_params, location)
        ! Read and validate simulation dates and calendar
        call validate_input_dates(cfg_params, start_datetime, end_datetime, calendar)  
        ! Load data every year or full load at the beginning.   
        load_yearly = cfg_params%get_param_logical('data.load_yearly', default=.false.)    

        ! Print header for simulation
        call print_header(location,start_datetime,end_datetime, load_yearly) 

        ! Verify if biogeochemistry is enabled
        is_bio_enabled = cfg_params%get_param_logical('biogeochemistry.enabled', default=.false.)
         if (is_bio_enabled) then
            is_sed_enabled = cfg_params%get_param_logical('biogeochemistry.sediments.enabled', default=.false.)
        end if

        ! Build vertical grids: water, sediment (if enabled), and full grid
        call build_grids(cfg_params, location%depth, is_bio_enabled, is_sed_enabled, wat_grid, sed_grid, full_grid, static_profs)

        ! Verify and initialise forcing data
        call PhysForc%set_error_mode(stop_on_error)

        call PhysForc%init(cfg_params, calendar, location, start_datetime, end_datetime, load_yearly, ok, msg)
        if (.not. ok) stop 'init_physics_forcing failed: '//trim(msg)
        ! Set calendar from forcing data
        if (calendar%kind == cal_unknown) then
            calendar%kind = PhysForc%get_sim_calendar()
        end if
        ! Calendar must be defined before initialising the simulation clock
        if (calendar%kind == cal_unknown) then
            stop 'Error: simulation calendar is unknown. '// &
                 'Set the calendar explicitly in the configuration or ensure the forcing NetCDF time variable has a valid CF calendar attribute.'
        end if

        ! Reference date is the start datetime for the simulation
        dt = cfg_params%get_param_int('time.dt', default=300, positive=.true.)        
        call init_clock(calendar,start_datetime, end_datetime, dt, sim_length_sec, n_steps, last_dt_length)

        call PhysForc%prepare(dt, ok=ok, errmsg=msg)
        if (.not. ok) stop 'prepare_physics_forcing failed: '//trim(msg)

        ! Load initial profiles to e.g. restart a simulation
        load_init_state = cfg_params%get_param_logical('load_initial_state.enabled', default=.false.)
        if (load_init_state) then            
            init_state_path = cfg_params%get_param_str('load_initial_state.filename', required=.true., trim_value=.true.)
            write(*,'(A)') 'Scanning initial state file: '//trim(init_state_path)

            call load_state_file(init_state_path, init_state, ok, msg)
            if (.not. ok) then
                stop 'Initial state loading failed: '//trim(msg)
            end if

            write(*,'(A,I0,A)') '  Found ', init_state%nvars, ' variables.'
        end if

        ! Initialise physics
        call init_physics(cfg_params, location, wat_grid, PE, load_init_state, init_state)

        ! Initialise biogeochemistry if it's the case
        if (is_bio_enabled) then
            call init_bio_fabm(cfg_params, location, wat_grid, sed_grid, full_grid,     &
                               start_datetime, end_datetime, calendar, load_yearly, dt, &
                               PE%PS, ForcSnp, PE%params%h0b, BE, EVT,                   &
                               load_init_state, init_state, static_profs)
            nsed = BE%nsed
        else 
            nsed = 0
            write(*,'(A)') 'Preparing a physics-only simulation (biogeochemistry is turned off).'
        end if

        if (is_bio_enabled) then
            all_vars = [PE%phys_vars, BE%env_int_vars, BE%int_vars, BE%btm_vars, BE%sfc_vars, &
                        BE%diag_hz_vars, BE%diag_int_vars, BE%conserved_vars,                 &
                        BE%tot_swiflux_vars, BE%dif_swiflux_vars, BE%bio_swiflux_vars]            
        else
            all_vars = [PE%phys_vars]
        end if
        ! Data needed for the output manager
        time_units = 'seconds since ' // trim(datetime_to_str(start_datetime)) ! Time units (CF-metadata convention)
        calname    = trim(calendar%name())

        ! Initialise output manager and output file
        call OM%init(cfg_params, full_grid, dt_s=dt, time_units=time_units, calendar_name=calname, &
                     vars=all_vars, loc=location, final_timestamp=datetime_to_str(end_datetime),   &
                     static_profiles=static_profs)

        is_main_initialized = .true.     
        
    end subroutine init_shelfseas

    !=============================================
    ! Main subroutine - Main loop over time is here
    !============================================
    subroutine run_shelfseas()
        implicit none

        integer(lk)        :: istep, model_time, dt_now
        integer            :: doy
        real(rk)           :: sec_of_day, doy_real
        integer            :: ierr
        character(len=512) :: errmsg
        logical            :: is_first_step = .true.  

        if (.not. is_main_initialized) error stop 'run_shelfseas: shelfseas not initialised.'

        model_time = 0_lk   ! Seconds since the start datetime of the simulation
       

        write(*,'(/,A,I0,A)') 'Starting simulation (', n_steps, ' steps)...'

        ! Main simulation loop
        do istep = 1_lk, n_steps
            if (istep > 1_lk) is_first_step = .false.
            ! Choose step length: dt for all but last step
            if (istep < n_steps) then
                dt_now = dt  ! Current time-step length, always the same except for the last step
            else
                dt_now = last_dt_length
            end if            

            ! Get current calendar time for FABM
            call simtime_to_datetime(calendar, start_datetime, model_time, current_datetime, doy)

            call PhysForc%tick(model_time, ok=ok, errmsg=msg)                               ! Time-manager to load yearly forcing data on time
            if (.not. ok) stop 'physics_forcing tick failed: '//trim(msg)
            call PhysForc%sample(model_time, ForcSnp, ok=ok, errmsg=msg)                    ! get forcing snapshot for the current model time
            if (.not. ok) stop 'physics_forcing sample failed: '//trim(msg)

            call solve_physics(PE, ForcSnp, dt_now, model_time, is_first_step, ierr, errmsg)

            if (ierr /= 0) then
                write(*,*) 'solve_physics failed: ', trim(errmsg)
                return
            end if

            if (is_bio_enabled) then
                sec_of_day = real(current_datetime%hour*3600 + current_datetime%minute*60 + current_datetime%second, rk)
                ! 0-based day-of-year + fractional day
                doy_real = real(doy - 1, rk) + sec_of_day / 86400._rk
                call integrate_bio_fabm(BE, PE%PS, ForcSnp, EVT, dt_now, istep, model_time, current_datetime, sec_of_day, doy_real)
                ! Update bioattenuation coefficients in Physics
                if (PE%params%apply_heat_bioshade) then
                    PE%atten_bio(:) = BE%BS%atten_coeff(BE%k_wat_btm:BE%k_wat_sfc)
                end if
            end if 
            ! Sample state for output
            call OM%step(dt_now)

            ! Progress bar
            call print_progress(istep, n_steps, model_time, sim_length_sec,  &
                                current_datetime, doy, start_datetime, end_datetime)
            ! Model time in seconds since start date
            model_time = model_time + dt_now
        end do
        write(*,'(A,/)') 'Simulation completed.'
    end subroutine run_shelfseas


    !======================
    ! Clean memory
    !======================
    subroutine end_shelfseas()
        if (.not. is_main_initialized) return
        call OM%close(sync_now=.true.)
        call end_physics(PE)
        if (is_bio_enabled) then
            call end_bio_fabm(BE)
        end if
        call PhysForc%clear()
        call clear_static_profiles(static_profs)
        call cfg_params%clear()
        ! Clear events
        call EVT%clear()
    end subroutine end_shelfseas 

end module shelfseas
