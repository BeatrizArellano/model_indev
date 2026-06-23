module event_manager
   use bio_types,        only: BioEnv
   use cf_time_utils,    only: seconds_between_datetimes
   use event_types,      only: Event, EVENT_TRACER_PULSE, EVENT_TRACER_REMOVAL, EVENT_TRAWLING
   use precision_types,  only: rk, lk
   use read_config_yaml, only: ConfigParams
   use str_utils,        only: to_lower
   use time_types,       only: DateTime, CFCalendar
   use time_utils,       only: datetime_from_string, datetime_to_str, parse_interval_to_seconds
   use tracer_pulse,     only: tracer_pulse_prepare, tracer_pulse_apply
   use tracer_removal,   only: tracer_removal_prepare, tracer_removal_apply

   implicit none
   private

   public :: EventManager

   character(len=16), parameter :: event_type_choices(3) = [character(len=16) :: &
                                                            EVENT_TRACER_PULSE, EVENT_TRACER_REMOVAL, EVENT_TRAWLING]

   type :: EventManager
      logical :: is_init          = .false.
      logical :: any_events       = .false.

      character(len=:), allocatable :: events_cfgfile

      type(ConfigParams) :: events_cfg
      type(Event),       allocatable :: events(:)
   contains
      procedure :: init                       => em_init
      procedure :: clear                      => em_clear
      procedure :: apply_tendencies           => em_apply_tendencies
      procedure :: apply_instantaneous        => em_apply_instantaneous
      procedure :: apply_direct_state_changes => em_apply_direct_state_changes
      procedure :: compute_next_dt            => em_compute_next_dt   
      !--- Private
      procedure, private :: read_events    => em_read_events
      procedure, private :: prepare_events => em_prepare_events
      procedure, private :: print_summary  => em_print_summary
   end type EventManager

contains

   subroutine em_init(self, event_cfg_file, BE, sim_startdate, sim_enddate, calendar)
      class(EventManager), intent(inout) :: self
      character(*),        intent(in)    :: event_cfg_file
      type(BioEnv),        intent(in)    :: BE
      type(DateTime),      intent(in)    :: sim_startdate, sim_enddate
      type(CFCalendar),    intent(in)    :: calendar

      call self%clear()

      if (len_trim(event_cfg_file) == 0 .or. trim(event_cfg_file) == 'off') then
         self%is_init = .true.
         return
      end if

      self%events_cfgfile = trim(event_cfg_file)

      call self%events_cfg%init()
      call self%events_cfg%load_yaml_content(self%events_cfgfile)

      call self%read_events(sim_startdate, sim_enddate, calendar)

      self%any_events = allocated(self%events) .and. size(self%events) > 0
      if (self%any_events) call self%prepare_events(BE)
      call self%print_summary()
      self%is_init = .true.
   end subroutine em_init  

   subroutine em_apply_tendencies(self, BE, model_time)
      class(EventManager), intent(in)    :: self
      type(BioEnv),        intent(inout) :: BE
      real(rk),            intent(in)    :: model_time

      integer :: ievent

      if (.not. self%any_events) return

      BE%tendency_int_evt = 0._rk

      do ievent = 1, size(self%events)         
         if (.not. self%events(ievent)%changes_tendency) cycle
         if (.not. self%events(ievent)%is_active(model_time)) cycle   ! Checks that the event is active during the current time-step
         if (self%events(ievent)%is_instantaneous) cycle

         select case (trim(self%events(ievent)%event_type))
            case (EVENT_TRACER_PULSE)
               call tracer_pulse_apply(self%events(ievent), BE)

            case default
               call fatal('event_manager:apply_tendencies', &
                          'Event is not implemented to be applied as tendency: '//trim(self%events(ievent)%event_type))
         end select
      end do
   end subroutine em_apply_tendencies

   ! Trigger events that produce instantaneous changes
   subroutine em_apply_instantaneous(self, BE, model_time)
      class(EventManager), intent(inout) :: self
      type(BioEnv),        intent(inout) :: BE
      real(rk),            intent(in)    :: model_time

      integer  :: ievent
      real(rk) :: eps_t

      if (.not. self%any_events) return

      eps_t = 1.0e-10_rk * max(1.0_rk, abs(model_time))

      do ievent = 1, size(self%events)

         if (.not. self%events(ievent)%is_instantaneous) cycle
         if (.not. self%events(ievent)%should_trigger(model_time)) cycle
         if (self%events(ievent)%was_triggered) cycle

         if (abs(model_time - self%events(ievent)%t_start) > eps_t) cycle

         !--------------------------------------------------------
         !  Apply event
         !--------------------------------------------------------
         select case (trim(self%events(ievent)%event_type))
         case (EVENT_TRAWLING)
            ! call trawling_apply(model_time, self%events(ievent), BE)
            call fatal('event_manager:apply_instantaneous', &
                        'Trawling event selected, but trawling_apply is not connected yet.')

         case default
            call fatal('event_manager:apply_instantaneous', &
               'Event is not implemented to be applied as instantaneous: '//trim(self%events(ievent)%event_type))
         end select

         self%events(ievent)%was_triggered = .true.

      end do
   end subroutine em_apply_instantaneous

   ! Direct state-change events are applied at the end of the substep,
   ! but their amount is integrated over the interval [model_time, model_time + dt].
   ! EventManager ensures substeps land on event start/end boundaries.
   subroutine em_apply_direct_state_changes(self, BE, model_time, dt)
      class(EventManager), intent(inout) :: self
      type(BioEnv),        intent(inout) :: BE
      real(rk),            intent(in)    :: model_time
      real(rk),            intent(in)    :: dt

      integer  :: ievent
      real(rk) :: t0, t1, t0_evt, t1_evt, dt_evt
      real(rk) :: eps_t

      if (.not. self%any_events) return
      if (dt <= 0._rk) return

      t0 = model_time
      t1 = model_time + dt
      eps_t = 1.0e-10_rk * max(1.0_rk, abs(t0), abs(t1), dt)

      do ievent = 1, size(self%events)

         if (.not. self%events(ievent)%changes_state_directly) cycle
         if (self%events(ievent)%is_instantaneous) cycle

         t0_evt = max(t0, self%events(ievent)%t_start)
         t1_evt = min(t1, self%events(ievent)%t_end)
         dt_evt = t1_evt - t0_evt

         if (dt_evt <= eps_t) cycle

         select case (trim(self%events(ievent)%event_type))

         case (EVENT_TRACER_REMOVAL)
            call tracer_removal_apply(self%events(ievent), BE, t0_evt, dt_evt)

         case default
            call fatal('event_manager:apply_state_changes', 'Event is not implemented to change state directly: '// &
                        trim(self%events(ievent)%event_type))
         end select

      end do
   end subroutine em_apply_direct_state_changes


   ! Compute the needed time-step to land exactly on the start/end of an event
   subroutine em_compute_next_dt(self, model_time, dt_left, dt_event_safe)
      class(EventManager), intent(in) :: self
      real(rk), intent(in)  :: model_time, dt_left
      real(rk), intent(out) :: dt_event_safe

      integer  :: i
      real(rk) :: t_now, t_next, t_limit, eps_t

      dt_event_safe = dt_left
      if (dt_left <= 0._rk) return
      if (.not. self%any_events) return

      t_now   = model_time
      t_limit = t_now + dt_left

      do i = 1, size(self%events)

         if (self%events(i)%is_instantaneous) then

            if (.not. self%events(i)%was_triggered) then
               t_next = self%events(i)%t_start
               eps_t = 1.0e-9_rk * max(1.0_rk, abs(t_now), abs(t_next), dt_left)

               if (t_next > t_now + eps_t .and. t_next <= t_limit + eps_t) then
                  dt_event_safe = min(dt_event_safe, t_next - t_now)
               end if
            end if

         else

            ! Start boundary
            t_next = self%events(i)%t_start
            eps_t = 1.0e-9_rk * max(1.0_rk, abs(t_now), abs(t_next), dt_left)
            if (t_next > t_now + eps_t .and. t_next <= t_limit + eps_t) then
               dt_event_safe = min(dt_event_safe, t_next - t_now)
            end if

            ! End boundary
            t_next = self%events(i)%t_end
            eps_t = 1.0e-9_rk * max(1.0_rk, abs(t_now), abs(t_next), dt_left)
            if (t_next > t_now + eps_t .and. t_next <= t_limit + eps_t) then
               dt_event_safe = min(dt_event_safe, t_next - t_now)
            end if

         end if

      end do

   end subroutine em_compute_next_dt

   subroutine em_clear(self)
      class(EventManager), intent(inout) :: self

      integer :: i

      if (allocated(self%events)) then
         do i = 1, size(self%events)
            call self%events(i)%clear()
         end do
         deallocate(self%events)
      end if

      if (allocated(self%events_cfgfile)) deallocate(self%events_cfgfile)

      call self%events_cfg%clear()

      self%is_init          = .false.
      self%any_events       = .false.
   end subroutine em_clear


   !==================================================================
   !           Internal
   !==================================================================
   subroutine em_read_events(self, sim_startdate, sim_enddate, calendar)
      class(EventManager), intent(inout) :: self
      type(DateTime),      intent(in)    :: sim_startdate, sim_enddate
      type(CFCalendar),    intent(in)    :: calendar

      character(len=:), allocatable :: event_names(:)
      character(len=:), allocatable :: dates(:)
      character(len=:), allocatable :: date_str
      character(len=:), allocatable :: basekey, dateskey, typekey
      character(len=:), allocatable :: event_type
      character(len=:), allocatable :: durationkey
      character(len=:), allocatable :: duration_str
      character(len=:), allocatable :: enabledkey

      logical     :: has_duration, enabled
      integer(lk) :: duration_seconds, sim_duration, event_start_sec
      integer     :: i, idate, ievent
      integer     :: total_events
      type(DateTime) :: event_datetime

      logical            :: ok
      character(len=256) :: errmsg


      if (.not. self%events_cfg%has_key('events')) then
         write(*,'(A)') &
            'WARNING: Event configuration file was loaded successfully, but no "events" section was found. Expected structure: events: <event_name>: ...'
         return
      end if

      sim_duration = seconds_between_datetimes(calendar, sim_enddate, sim_startdate)
      event_names = self%events_cfg%get_child_keys('events')

      !--------------------------------------------------
      ! Count total internal events depending on dates
      !--------------------------------------------------
      total_events = 0

      do i = 1, size(event_names)
         basekey  = 'events.'//trim(event_names(i))
         dateskey = basekey//'.dates'
         enabledkey = basekey // '.enabled'

         if (self%events_cfg%has_key(enabledkey)) then
            enabled = self%events_cfg%get_param_logical(enabledkey, required=.true.)
         else
            enabled = .true.
         end if
         if (.not. enabled) cycle

         if (self%events_cfg%is_string_list(dateskey)) then
            dates = self%events_cfg%get_param_str_list(dateskey)
            total_events = total_events + size(dates)

         else if (self%events_cfg%is_string(dateskey)) then
            total_events = total_events + 1

         else
            call fatal('event_manager:read_events', &
               'Event "'//trim(event_names(i))//'" must define "dates" as a string or string list.')
         end if
      end do

      allocate(self%events(total_events))

      if (total_events == 0) then
         write(*,'(A)') &
            'WARNING: No events were found in the "events" section.'
         return
      end if

      !--------------------------------
      ! Populate expanded event list
      !-------------------------------
      ievent = 0

      do i = 1, size(event_names)
         basekey = 'events.'//trim(event_names(i))
         typekey = basekey//'.type'
         dateskey = basekey//'.dates'
         enabledkey = basekey // '.enabled'

         if (self%events_cfg%has_key(enabledkey)) then
            enabled = self%events_cfg%get_param_logical(enabledkey, required=.true.)
         else
            enabled = .true.
         end if
         if (.not. enabled) cycle

         event_type = self%events_cfg%get_param_str(typekey, required=.true., choices=event_type_choices, trim_value=.true., match_case=.false.)

         if (self%events_cfg%is_string_list(dateskey)) then
            dates = self%events_cfg%get_param_str_list(dateskey)
         else
            date_str = self%events_cfg%get_param_str(dateskey, required=.true., trim_value=.true.)
            if (allocated(dates)) deallocate(dates)
            allocate(character(len=len_trim(date_str)) :: dates(1))
            dates(1) = trim(date_str)
         end if

         durationkey = basekey//'.duration'

         has_duration = self%events_cfg%has_key(durationkey)
         if (has_duration) then
            duration_str = self%events_cfg%get_param_str(durationkey, required=.true., trim_value=.true.)
            call parse_interval_to_seconds(duration_str, duration_seconds)   ! e.g., "1h" -> 3600
         else
            duration_seconds = 0_lk
         end if

         do idate = 1, size(dates)
            ievent = ievent + 1

            self%events(ievent)%name             = trim(event_names(i))
            self%events(ievent)%event_type       = to_lower(trim(event_type))
            self%events(ievent)%date_string      = trim(dates(idate))
            self%events(ievent)%is_instantaneous = .not. has_duration
            self%events(ievent)%definition_idx   = i
            self%events(ievent)%occurrence_idx   = idate

            event_datetime = datetime_from_string(trim(dates(idate)), ok, errmsg, calendar%kind)
            if (.not. ok) then
               call fatal('event_manager:read_events', &
                  'Invalid event date "'//trim(dates(idate))//'": '//trim(errmsg))
            end if
            event_start_sec = seconds_between_datetimes(calendar, event_datetime, sim_startdate)

            if (event_start_sec < 0_lk) then
               call fatal('event_manager:read_events', &
                  'Event "'//trim(event_names(i))//'" occurs before the simulation start: '//trim(dates(idate)))
            end if

            if (event_start_sec >= sim_duration) then
               call fatal('event_manager:read_events', &
                  'Event "'//trim(event_names(i))//'" occurs at or after the simulation end: '//trim(dates(idate)))
            end if

            self%events(ievent)%t_start = real(event_start_sec, rk)

            if (has_duration) then
               if (duration_seconds <= 0_lk) then
                  call fatal('event_manager:read_events', &
                     'Event "'//trim(event_names(i))//'" has invalid duration; use positive forms like 900s, 15m, 1h, 1d.')
               end if

               self%events(ievent)%duration = real(duration_seconds, rk)
               self%events(ievent)%t_end = self%events(ievent)%t_start + real(duration_seconds, rk)

               if (self%events(ievent)%t_end > real(sim_duration, rk)) then
                  write(*,'(A)') 'WARNING: Event "'//trim(event_names(i))//'" extends beyond the end of simulation.'
               end if

            else
               self%events(ievent)%duration = 0._rk
               self%events(ievent)%t_end = self%events(ievent)%t_start
            end if
         end do
      end do

   end subroutine em_read_events

   subroutine em_prepare_events(self, BE)
      class(EventManager), intent(inout) :: self
      type(BioEnv),        intent(in)    :: BE

      integer :: i

      if (.not. allocated(self%events)) return

      do i = 1, size(self%events)

         select case (trim(self%events(i)%event_type))

            case (EVENT_TRACER_PULSE)
               call tracer_pulse_prepare(self%events(i), self%events_cfg, BE)

            case (EVENT_TRACER_REMOVAL)
               call tracer_removal_prepare(self%events(i), self%events_cfg, BE)

            case (EVENT_TRAWLING)
               !call trawling_prepare_event(self%events(i), self%events_cfg)

            case default
               call fatal('event_manager:prepare_events', &
                  'Unknown event type: '//trim(self%events(i)%event_type))

         end select

      end do
   end subroutine em_prepare_events


   subroutine em_print_summary(self)
      class(EventManager), intent(in) :: self

      integer :: i, j, ndefs, nocc
      logical :: first_date

      if (.not. allocated(self%events)) then
         write(*,'(A)') 'A configuration file for external events was found, but no active events were detected.'
         return
      end if

      if (size(self%events) == 0) then
         write(*,'(A)') 'A configuration file for external events was found, but no active events were detected.'
         return
      end if

      ndefs = 0
      do i = 1, size(self%events)
         if (self%events(i)%occurrence_idx == 1) ndefs = ndefs + 1
      end do

      nocc = size(self%events)

      write(*,'(A)')    'External events:'
      write(*,'(A,I0, A, I0)') ' Event definitions: ', ndefs, ', Total occurrences: ', nocc
      write(*,'(A)')    ''

      do i = 1, size(self%events)

         if (self%events(i)%occurrence_idx /= 1) cycle

         write(*,'(A,A,A,A,A)') ' -', trim(self%events(i)%name), &
            ' [', trim(self%events(i)%event_type), ']'

         if (self%events(i)%is_instantaneous) then
            write(*,'(A)') '    mode: instantaneous'
         else
            write(*,'(A,F0.2,A)') '    mode: interval, duration: ', &
               self%events(i)%duration, ' s'
         end if

         write(*,'(A)', advance='no') '    occurrences: '

         first_date = .true.

         do j = 1, size(self%events)
            if (self%events(j)%definition_idx /= self%events(i)%definition_idx) cycle

            if (.not. first_date) write(*,'(A)', advance='no') ', '
            write(*,'(A)', advance='no') trim(self%events(j)%date_string)

            first_date = .false.
         end do

         write(*,*)
         write(*,*)

      end do
   end subroutine em_print_summary


   subroutine fatal(where, msg)
      character(len=*), intent(in) :: where, msg
      write(*,'(A,1X,A)') '[FATAL:'//trim(where)//']', trim(msg)
      stop 1
   end subroutine fatal

end module event_manager