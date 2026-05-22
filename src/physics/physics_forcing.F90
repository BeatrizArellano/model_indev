! src/physics/physics_forcing.F90
!
! Adapter between the generic input-data layer and the
! physics core. This module owns the physics forcing DataManager, reads the
! physics forcing YAML section, builds generic DataSpec(:), and returns a
! ForcingSnapshot for each model step.
!
module physics_forcing
   use data_manager,     only: DataManager
   use data_types,       only: DataLoaderCfg, DataSpec, &
                               DATA_INPUT_FILE, DATA_INPUT_CONSTANT, DATA_INPUT_OFF, &
                               DATA_TIME_ABSOLUTE, DATA_TIME_REPEAT_YEAR
   use geo_utils,        only: LocationInfo   
   use precision_types,  only: rk, lk
   use read_config_yaml, only: ConfigParams
   use str_utils,        only: to_lower, inttostr
   use time_types,       only: DateTime, CFCalendar, calendar_compatible, cal_gregorian, cal_unknown
   use time_utils,       only: is_leap_gregorian

   implicit none
   private

   public :: ForcingSnapshot
   public :: PhysicsForcing

   integer, parameter :: N_PHYS_FORCING = 10

   character(len=24), parameter :: phys_forcing_names(N_PHYS_FORCING) = [ character(len=24) :: &
                                                                        'surf_air_temp',       &
                                                                        'sl_pressure',         &
                                                                        'relative_humidity',   &
                                                                        'shortwave_radiation', &
                                                                        'longwave_radiation',  &
                                                                        'wind_u10',            &
                                                                        'wind_v10',            &
                                                                        'precipitation',       &
                                                                        'evaporation',         &
                                                                        'runoff'               &
                                                                    ]


   ! Forcing snapshot.
   ! This stores forcing fields loaded from file or constants.
   type :: ForcingSnapshot
      real(rk) :: air_temp  = 0.0_rk
      real(rk) :: slp       = 0.0_rk
      real(rk) :: rel_hum   = 0.0_rk
      real(rk) :: short_rad = 0.0_rk
      real(rk) :: long_rad  = 0.0_rk

      real(rk) :: wind_u10  = 0.0_rk   ! zonal wind component at 10 m [m s-1]
      real(rk) :: wind_v10  = 0.0_rk   ! meridional wind component at 10 m [m s-1]

      ! Freshwater fluxes [m s-1]. 
      real(rk) :: precip    = 0.0_rk
      real(rk) :: evap      = 0.0_rk
      real(rk) :: runoff    = 0.0_rk
   end type ForcingSnapshot


   type :: PhysicsForcing
      type(DataManager) :: dm

      logical :: is_init       = .false.
      logical :: is_prepared   = .false.
      logical :: stop_on_error = .true.

   contains
      procedure :: init
      procedure :: prepare
      procedure :: tick
      procedure :: sample
      procedure :: clear
      procedure :: set_error_mode
      procedure :: get_sim_calendar
   end type PhysicsForcing

contains


   subroutine init(self, params, calendar_cfg, location, start_datetime, end_datetime, ok, errmsg)
      class(PhysicsForcing), intent(inout) :: self
      type(ConfigParams),    intent(in)    :: params
      type(CFCalendar),      intent(in)    :: calendar_cfg
      type(LocationInfo),    intent(in)    :: location
      type(DateTime),        intent(in)    :: start_datetime, end_datetime
      logical,               intent(out)   :: ok
      character(*),          intent(out)   :: errmsg

      type(DataLoaderCfg)        :: cfg
      type(DataSpec), allocatable :: specs(:)
      logical :: lok
      character(len=512) :: lmsg
      logical :: has_file_forcing
      integer :: i

      ok = .false.
      errmsg = ''

      call self%clear()
      call self%dm%set_error_mode(self%stop_on_error)

      !-----------------------------
      ! Read forcing parameters. 
      !-----------------------------
      call read_physics_forcing_config(params, calendar_cfg, start_datetime, end_datetime, specs, cfg, lok, lmsg)
      if (.not. lok) then
         errmsg = trim(lmsg)
         call stop_fatal('init/read_physics_forcing_config', errmsg, self%stop_on_error)
         return
      end if

      !------------------------------------------------------------
      ! If calendar=0, the calendar must be derived from file data.
      ! This is impossible when all active forcing variables are
      ! constants or off, so fail before scanning/loading any file.
      !------------------------------------------------------------
      has_file_forcing = .false.
      do i = 1, size(specs)
         if (specs(i)%input_type == DATA_INPUT_FILE) then
            has_file_forcing = .true.
            exit
         end if
      end do

      if (calendar_cfg%kind == cal_unknown .and. .not. has_file_forcing) then
         errmsg = 'time.calendar=0 requests deriving the calendar from forcing data, ' // &
                  'but no active physics forcing variable is read from file. ' // &
                  'Set time.calendar explicitly, or set at least one forcing variable to mode=file.'
         call stop_fatal('init/calendar', errmsg, self%stop_on_error)
         return
      end if

      !-----------------------------
      ! Scan, check and prepare forcing data. 
      !-----------------------------
      call self%dm%init(specs, cfg, calendar_cfg, location, start_datetime, end_datetime, lok, lmsg)
      if (.not. lok) then
         errmsg = trim(lmsg)
         call stop_fatal('init/DataManager%init', errmsg, self%stop_on_error)
         return
      end if

      if (self%dm%sim_cal%kind == cal_unknown) then
         errmsg = 'Physics forcing calendar is still unknown after DataManager initialization.'
         call stop_fatal('init/calendar', errmsg, self%stop_on_error)
         return
      end if

      ! In case a year is repeated 
      if (cfg%repeat_enabled) then
         write(*,'(a,i0,a)') &
            'Physics forcing: repeat-year mode enabled. Forcing data from year ', &
            cfg%repeat_year, ' will be reused for all simulation years.'
         if (calendar_compatible(self%dm%sim_cal%kind, cal_gregorian) .and. &
             .not. is_leap_gregorian(cfg%repeat_year)) then

            write(*,'(a)') &
               ' Warning: Repeating a non-leap Gregorian year. Feb 29 forcing will use Feb 28 values.'
         end if

      end if

      self%is_init     = .true.
      self%is_prepared = .false.

      ok = .true.
   end subroutine init


   subroutine prepare(self, dt_main, ok, errmsg)
      class(PhysicsForcing), intent(inout) :: self
      integer(lk),           intent(in)    :: dt_main
      logical,               intent(out)   :: ok
      character(*),          intent(out)   :: errmsg

      ok = .false.
      errmsg = ''

      if (.not. self%is_init) then
         errmsg = 'PhysicsForcing not initialized'
         call stop_fatal('prepare', errmsg, self%stop_on_error)
         return
      end if

      !-------------------------
      ! Load data
      !----------------------------
      call self%dm%prepare(dt_main, ok, errmsg)
      if (.not. ok) then
         call stop_fatal('prepare/DataManager%prepare', errmsg, self%stop_on_error)
         return
      end if

      self%is_prepared = .true.
   end subroutine prepare

   !-------------------------------------------------
   ! Updating forcing data to be ready when requested
   !-------------------------------------------------
   subroutine tick(self, model_time, ok, errmsg)
      class(PhysicsForcing), intent(inout) :: self
      integer(lk),           intent(in)    :: model_time
      logical, optional,     intent(out)   :: ok
      character(*), optional,intent(out)   :: errmsg

      logical :: lok
      character(len=512) :: lmsg

      if (present(ok))     ok = .true.
      if (present(errmsg)) errmsg = ''

      if (.not. self%is_prepared) return

      call self%dm%tick(model_time, lok, lmsg)
      if (.not. lok) then
         if (present(ok))     ok = .false.
         if (present(errmsg)) errmsg = trim(lmsg)
         call stop_fatal('tick/DataManager%tick', lmsg, self%stop_on_error)
      end if
   end subroutine tick

   !-------------------------------------------------------
   ! Get forcing snapshot
   !------------------------------------------------------
   subroutine sample(self, model_time, fs, ok, errmsg)
      class(PhysicsForcing), intent(inout) :: self
      integer(lk),           intent(in)    :: model_time
      type(ForcingSnapshot), intent(out)   :: fs
      logical, optional,     intent(out)   :: ok
      character(*), optional,intent(out)   :: errmsg

      logical :: lok
      character(len=512) :: lmsg

      if (present(ok))     ok = .true.
      if (present(errmsg)) errmsg = ''

      if (.not. self%is_prepared) then
         call fail_sample('sample', 'PhysicsForcing not prepared', self%stop_on_error, ok, errmsg)
         return
      end if

      fs%air_temp = self%dm%value('surf_air_temp', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%slp = self%dm%value('sl_pressure', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%rel_hum = self%dm%value('relative_humidity', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%short_rad = self%dm%value('shortwave_radiation', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%long_rad = self%dm%value('longwave_radiation', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%wind_u10 = self%dm%value('wind_u10', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%wind_v10 = self%dm%value('wind_v10', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%precip = self%dm%value('precipitation', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%evap = self%dm%value('evaporation', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

      fs%runoff = self%dm%value('runoff', model_time, lok, lmsg)
      if (.not. lok) then; call fail_sample('sample/DataManager%value', lmsg, self%stop_on_error, ok, errmsg); return; end if

   end subroutine sample


   subroutine clear(self)
      class(PhysicsForcing), intent(inout) :: self

      call self%dm%clear()
      self%is_init     = .false.
      self%is_prepared = .false.
   end subroutine clear


   subroutine set_error_mode(self, on)
      class(PhysicsForcing), intent(inout) :: self
      logical,               intent(in)    :: on

      self%stop_on_error = on
      call self%dm%set_error_mode(on)
   end subroutine set_error_mode


   integer function get_sim_calendar(self) result(cal_kind)
      class(PhysicsForcing), intent(in) :: self 

       cal_kind = self%dm%get_sim_calendar()
   end function get_sim_calendar
 

    !--------------------------------------
    !  Internal helpers
    !---------------------------------------
    subroutine fail_sample(where, message, stop_on_error, ok, errmsg)
        character(*), intent(in)            :: where
        character(*), intent(in)            :: message
        logical,      intent(in)            :: stop_on_error
        logical,      intent(out), optional :: ok
        character(*), intent(out), optional :: errmsg

        if (present(ok))     ok = .false.
        if (present(errmsg)) errmsg = trim(message)

        call stop_fatal(where, message, stop_on_error)
    end subroutine fail_sample

   subroutine read_physics_forcing_config(params, calendar, start_datetime, end_datetime, specs, cfg, ok, errmsg)
      type(ConfigParams),        intent(in)  :: params
      type(CFCalendar),          intent(in)  :: calendar
      type(DateTime),            intent(in)  :: start_datetime, end_datetime
      type(DataSpec), allocatable, intent(out) :: specs(:)
      type(DataLoaderCfg),       intent(out) :: cfg
      logical,                   intent(out) :: ok
      character(*),              intent(out) :: errmsg

      character(:), allocatable :: global_file
      character(:), allocatable :: sal_mode
      character(len=8), dimension(2) :: sal_choices
      integer :: i, ip, ie

      ok = .false.
      errmsg = ''

      sal_choices  = ['constant', 'compute ']

      allocate(specs(N_PHYS_FORCING))

      cfg%cfg_calendar = calendar%kind
      cfg%load_yearly  = params%get_param_logical('forcing.load_yearly', default=.false., strict=.false.)
      cfg%time_mode    = DATA_TIME_ABSOLUTE

      if (.not. params%is_disabled('forcing.repeat_forcing_year')) then
         cfg%repeat_year    = params%get_param_int('forcing.repeat_forcing_year')
         cfg%repeat_enabled = .true.
         cfg%time_mode      = DATA_TIME_REPEAT_YEAR
      else
         cfg%repeat_enabled = .false.
         cfg%repeat_year    = -huge(1)
      end if

      global_file = normalise_optional_string(params%get_param_str('forcing.filename', default='', empty_ok=.true.))

      do i = 1, N_PHYS_FORCING
         call read_one_forcing_var(params, trim(phys_forcing_names(i)), global_file, specs(i), ok, errmsg)
         if (.not. ok) return
      end do

      ! Freshwater forcing is loaded only if salinity is prognostic/computed.
      sal_mode = to_lower(trim(params%get_param_str('physics.variables.salinity.mode', &
                          choices=sal_choices, trim_value=.true., match_case=.false., default='constant')))

      if (sal_mode == 'constant') then
         call force_spec_off(specs, 'precipitation')
         call force_spec_off(specs, 'evaporation')
         call force_spec_off(specs, 'runoff')
      else if (sal_mode == 'compute') then
         ip = find_spec_index(specs, 'precipitation')
         ie = find_spec_index(specs, 'evaporation')

         if (ip <= 0 .or. specs(ip)%input_type == DATA_INPUT_OFF) then
            errmsg = 'salinity mode=compute but forcing.precipitation is off or missing. '// &
                     'Set forcing.precipitation.mode to file or constant, or use salinity.mode=constant.'
            return
         end if

         if (ie <= 0 .or. specs(ie)%input_type == DATA_INPUT_OFF) then
            errmsg = 'salinity mode=compute but forcing.evaporation is off or missing. '// &
                     'Set forcing.evaporation.mode to file or constant, or use salinity.mode=constant.'
            return
         end if
      end if

      ok = .true.
   end subroutine read_physics_forcing_config


   subroutine read_one_forcing_var(params, id, global_file, spec, ok, errmsg)
      type(ConfigParams), intent(in)    :: params
      character(*),       intent(in)    :: id
      character(*),       intent(in)    :: global_file
      type(DataSpec),     intent(inout) :: spec
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      character(:), allocatable :: mode, source_name, file_name, time_name
      character(len=8), dimension(5) :: mode_choices
      logical :: required_var

      ok = .false.
      errmsg = ''

      mode_choices = ['file    ', 'constant', 'compute ', 'off     ', 'false   ']

      spec%name = trim(id)
      spec%units = ''

      required_var = is_required_physics_forcing(id)

      if (required_var) then
         mode = params%get_param_str('forcing.'//trim(id)//'.mode', required=.true., &
                                     choices=mode_choices, trim_value=.true., match_case=.false.)
      else
         mode = params%get_param_str('forcing.'//trim(id)//'.mode', default='off', &
                                     choices=mode_choices, trim_value=.true., match_case=.false.)
      end if

      mode = to_lower(trim(mode))

      select case (mode)

      case ('file')
         source_name = params%get_param_str('forcing.'//trim(id)//'.name', required=.true.)
         file_name   = params%get_param_str('forcing.'//trim(id)//'.filename', default='', empty_ok=.true.)
         time_name   = params%get_param_str('forcing.'//trim(id)//'.time_name', default='time')

         file_name = normalise_optional_string(file_name)

         if (len_trim(file_name) == 0) then
            if (len_trim(global_file) == 0) then
               errmsg = 'forcing.'//trim(id)//' has mode=file but no per-variable filename and no forcing.filename.'
               return
            end if
            file_name = trim(global_file)
         end if

         spec%input_type = DATA_INPUT_FILE
         spec%source_var = trim(source_name)
         spec%path       = trim(file_name)
         spec%time_var   = trim(time_name)

      case ('constant')
         spec%input_type  = DATA_INPUT_CONSTANT
         spec%const_value = params%get_param_num('forcing.'//trim(id)//'.constant', finite=.true.)
         spec%source_var  = ''
         spec%path        = ''
         spec%time_var    = ''

      case ('compute')
         errmsg = 'forcing.'//trim(id)//'.mode=compute is not implemented in physics_forcing yet.'
         return

      case ('off', 'false')
         if (required_var) then
            errmsg = 'Required physics forcing variable '//trim(id)//' cannot be off.'
            return
         end if

         spec%input_type = DATA_INPUT_OFF
         spec%source_var = ''
         spec%path       = ''
         spec%time_var   = ''

      case default
         errmsg = 'Invalid forcing mode for '//trim(id)//': '//trim(mode)
         return
      end select

      ok = .true.
   end subroutine read_one_forcing_var


   pure logical function is_required_physics_forcing(id) result(required)
      character(*), intent(in) :: id

      select case (trim(id))
      case ('surf_air_temp', 'sl_pressure', 'relative_humidity', &
            'shortwave_radiation', 'longwave_radiation', &
            'wind_u10', 'wind_v10')
         required = .true.
      case default
         required = .false.
      end select
   end function is_required_physics_forcing


   subroutine force_spec_off(specs, name)
      type(DataSpec), intent(inout) :: specs(:)
      character(*),   intent(in)    :: name

      integer :: i

      i = find_spec_index(specs, name)
      if (i <= 0) return

      specs(i)%input_type = DATA_INPUT_OFF
      specs(i)%source_var = ''
      specs(i)%path       = ''
      specs(i)%time_var   = ''
      specs(i)%const_value = 0.0_rk
   end subroutine force_spec_off


   integer function find_spec_index(specs, name) result(idx)
      type(DataSpec), intent(in) :: specs(:)
      character(*),   intent(in) :: name

      integer :: i

      idx = 0
      do i = 1, size(specs)
         if (.not. allocated(specs(i)%name)) cycle
         if (trim(specs(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do
   end function find_spec_index


   function normalise_optional_string(s) result(out)
      character(*), intent(in) :: s
      character(:), allocatable :: out

      character(:), allocatable :: low

      out = trim(adjustl(s))
      low = to_lower(out)

      if (low == 'off' .or. low == 'false' .or. low == 'null' .or. low == 'none' .or. low == '~') then
         out = ''
      end if
   end function normalise_optional_string


   subroutine stop_fatal(where, msg, stop_on_error)
      character(*), intent(in) :: where, msg
      logical,      intent(in) :: stop_on_error

      if (.not. stop_on_error) return

      write(*,*) 'FATAL PhysicsForcing['//trim(where)//']: ', trim(msg)
      stop 1
   end subroutine stop_fatal

end module physics_forcing