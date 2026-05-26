! src/io/data_manager.F90
module data_manager
   use cf_time_utils,   only: seconds_between_datetimes, seconds_since_datetime_file
   use data_loader,     only: DataLoaderState, scan_and_init_data, load_input_data, &
                              print_data_summary, find_data_index, init_series_cursor, &
                              value_at_step
   use data_types,      only: DataLoaderCfg, DataSpec, InputData, DataVarSeries, &
                              DATA_TIME_ABSOLUTE, DATA_TIME_REPEAT_YEAR, DATA_INPUT_FILE  
   use geo_utils,       only: LocationInfo   
   use precision_types, only: rk, lk
   use sim_clocks,      only: year_to_simyear, simtime_to_datetime
   use str_utils,       only: inttostr
   use time_types,      only: DateTime, CFCalendar, calendar_compatible, cal_gregorian, cal_unknown
   use time_utils,      only: is_leap_gregorian

   implicit none
   private

   public :: DataManager
   public :: compute_switch_time

   type :: DataManager
      type(DataLoaderState) :: state
      type(DateTime)        :: sim_start
      type(CFCalendar)      :: sim_cal

      type(InputData), allocatable :: data_curr
      type(InputData), allocatable :: data_next

      integer(lk) :: dt_main     = 0_lk
      integer(lk) :: switch_time = huge(1_lk)
      integer     :: y_active    = 0

      logical :: is_init       = .false.
      logical :: have_curr     = .false.
      logical :: have_next     = .false.
      logical :: stop_on_error = .true.

   contains
      procedure :: init
      procedure :: prepare
      procedure :: tick
      procedure :: value
      procedure :: sampling_time
      procedure :: get_index
      procedure :: clear      
      procedure :: set_error_mode
      procedure :: get_sim_calendar
   end type DataManager

contains

   subroutine init(self, specs, cfg, calendar_cfg, location, start_datetime, end_datetime, ok, errmsg)
      class(DataManager), intent(inout) :: self
      type(DataSpec),     intent(in)    :: specs(:)
      type(DataLoaderCfg),intent(in)    :: cfg
      type(CFCalendar),   intent(in)    :: calendar_cfg
      type(LocationInfo), intent(in)    :: location
      type(DateTime),     intent(in)    :: start_datetime, end_datetime
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      ok = .false.
      errmsg = ''

      write(*,'(A)') 'Scanning input data...'

      call self%clear()

      call scan_and_init_data(specs, cfg, calendar_cfg, location, start_datetime, end_datetime, &
                              self%state, ok, errmsg)

      if (.not. ok) then
         call stop_fatal('init/scan_and_init_data', errmsg, self%stop_on_error)
         return
      end if

      call print_data_summary(self%state)

      self%sim_cal%kind = self%state%sim_cal%kind
      self%sim_start    = start_datetime

      self%is_init   = .true.
      self%have_curr = .false.
      self%have_next = .false.

      ok = .true.
   end subroutine init


   subroutine prepare(self, dt_main, ok, errmsg)
      class(DataManager), intent(inout) :: self
      integer(lk),        intent(in)    :: dt_main
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      ok = .false.
      errmsg = ''

      if (.not. self%is_init) then
         errmsg = 'DataManager not initialized'
         call stop_fatal('prepare', errmsg, self%stop_on_error)
         return
      end if

      self%dt_main = dt_main

      if (self%state%cfg%load_yearly) then
         call prepare_yearly(self, ok, errmsg)
      else
         call prepare_full(self, ok, errmsg)
      end if

      if (.not. ok) call stop_fatal('prepare', errmsg, self%stop_on_error)
   end subroutine prepare


   subroutine prepare_yearly(self, ok, errmsg)
      class(DataManager), intent(inout) :: self
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      integer :: k
      logical :: lok
      character(len=512) :: lmsg

      ok = .false.
      errmsg = ''

      self%y_active    = self%sim_start%year
      self%switch_time = compute_switch_time(self%sim_cal, self%sim_start, self%y_active + 1)

      k = year_to_simyear(self%y_active, self%state%sim_y_start, self%state%sim_y_end)

      if (k <= 0) then
         errmsg = 'prepare_yearly: start year outside simulation range'
         return
      end if

      if (.not. allocated(self%data_curr)) allocate(self%data_curr)

      call load_input_data(self%state, k, self%data_curr, lok, lmsg)

      if (.not. lok) then
         errmsg = '(Y='//trim(adjustl(inttostr(self%sim_start%year)))//') failed: '//trim(lmsg)
         return
      end if

      self%have_curr = .true.
      self%have_next = .false.

      ok = .true.
   end subroutine prepare_yearly


   subroutine prepare_full(self, ok, errmsg)
      class(DataManager), intent(inout) :: self
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      logical :: lok
      character(len=512) :: lmsg

      ok = .false.
      errmsg = ''

      if (.not. allocated(self%data_curr)) allocate(self%data_curr)

      call load_full_data(self%state, self%data_curr, lok, lmsg)

      if (.not. lok) then
         errmsg = 'prepare_full: full input data load failed: '//trim(lmsg)
         return
      end if

      self%y_active    = self%sim_start%year
      self%switch_time = huge(1_lk)
      self%have_curr   = .true.
      self%have_next   = .false.

      ok = .true.
   end subroutine prepare_full


   subroutine tick(self, model_time, ok, errmsg)
      class(DataManager), intent(inout) :: self
      integer(lk),        intent(in)    :: model_time
      logical, optional,  intent(out)   :: ok
      character(*), optional, intent(out) :: errmsg

      integer :: k
      logical :: lok
      character(len=512) :: lmsg

      if (present(ok))     ok = .true.
      if (present(errmsg)) errmsg = ''

      if (.not. self%is_init)   return
      if (.not. self%have_curr) return

      ! Full-memory mode has no yearly switching.
      if (.not. self%state%cfg%load_yearly) return

      ! Preload next year one main time step before the boundary.
      if ((.not. self%have_next) .and. model_time >= self%switch_time - self%dt_main) then
         k = year_to_simyear(self%y_active + 1, self%state%sim_y_start, self%state%sim_y_end)

         if (k > 0) then
            if (.not. allocated(self%data_next)) allocate(self%data_next)

            call copy_loaded_repeat_series(self%state, self%data_curr, self%data_next)
            call load_input_data(self%state, k, self%data_next, lok, lmsg, skip_loaded_repeats=.true.)

            if (.not. lok) then
               if (present(ok))     ok = .false.
               if (present(errmsg)) errmsg = 'tick: preload next year failed: '//trim(lmsg)
               call stop_fatal('tick/preload', lmsg, self%stop_on_error)
               return
            end if

            self%have_next = .true.
         end if
      end if

      ! Promote next year at or after the year boundary.
      if (model_time >= self%switch_time) then
         if (self%have_next) then
            if (allocated(self%data_curr)) deallocate(self%data_curr)
            call move_alloc(from=self%data_next, to=self%data_curr)

            self%have_curr = .true.
            self%have_next = .false.
         else
            k = year_to_simyear(self%y_active + 1, self%state%sim_y_start, self%state%sim_y_end)

            if (k > 0) then
               if (.not. allocated(self%data_curr)) allocate(self%data_curr)

               call load_input_data(self%state, k, self%data_curr, lok, lmsg, skip_loaded_repeats=.true.)

               if (.not. lok) then
                  if (present(ok))     ok = .false.
                  if (present(errmsg)) errmsg = 'tick: just-in-time year load failed: '//trim(lmsg)
                  call stop_fatal('tick/just_in_time', lmsg, self%stop_on_error)
                  return
               end if

               self%have_curr = .true.
            end if
         end if

         self%y_active    = self%y_active + 1
         self%switch_time = compute_switch_time(self%sim_cal, self%sim_start, self%y_active + 1)
      end if
   end subroutine tick


   real(rk) function value(self, name, model_time, ok, errmsg) result(v)
      class(DataManager), intent(inout) :: self
      character(*),       intent(in)    :: name
      integer(lk),        intent(in)    :: model_time
      logical, optional,  intent(out)   :: ok
      character(*), optional, intent(out) :: errmsg

      integer :: idx

      v = 0.0_rk

      if (present(ok))     ok = .true.
      if (present(errmsg)) errmsg = ''

      if (.not. self%have_curr) then
         if (present(ok))     ok = .false.
         if (present(errmsg)) errmsg = 'DataManager%value: data are not loaded'
         call stop_fatal('value', 'data are not loaded', self%stop_on_error)
         return
      end if

      idx = find_data_index(self%data_curr, name)

      if (idx <= 0) then
         if (present(ok))     ok = .false.
         if (present(errmsg)) errmsg = 'DataManager%value: missing variable '//trim(name)
         call stop_fatal('value', 'missing variable '//trim(name), self%stop_on_error)
         return
      end if

      v = value_at_step(self%data_curr%vars(idx), self%sampling_time(self%data_curr%vars(idx), model_time))
   end function value

   integer(lk) function sampling_time(self, series, model_time) result(t_sample)
      class(DataManager), intent(in) :: self
      type(DataVarSeries), intent(in) :: series
      integer(lk),        intent(in) :: model_time

      type(DateTime) :: dt_model
      integer :: doy

      if (series%is_const) then
         t_sample = model_time
         return
      end if

      select case (series%time_mode)

      case (DATA_TIME_ABSOLUTE)

         ! Convert model-relative seconds into this series' source-file time axis.
         t_sample = model_time + series%sim_offset

      case (DATA_TIME_REPEAT_YEAR)

         call simtime_to_datetime(self%sim_cal, self%sim_start, model_time, dt_model, doy)

         if (calendar_compatible(series%cal%kind, cal_gregorian) .and. &
            dt_model%month == 2 .and. dt_model%day == 29 .and. &
            .not. is_leap_gregorian(series%repeat_year)) then

            t_sample = int(nint(seconds_since_datetime_file(series%cal, series%u, &
                                                            series%repeat_year, &
                                                            2, 28, &
                                                            dt_model%hour, dt_model%minute, dt_model%second)), lk)
         else
            t_sample = int(nint(seconds_since_datetime_file(series%cal, series%u, &
                                                            series%repeat_year, &
                                                            dt_model%month, dt_model%day, &
                                                            dt_model%hour, dt_model%minute, dt_model%second)), lk)
         end if

      case default

         ! Safe fallback: behave as transient/absolute forcing.
         t_sample = model_time + series%sim_offset

      end select
   end function sampling_time


   integer function get_index(self, name) result(idx)
      class(DataManager), intent(in) :: self
      character(*),       intent(in) :: name

      idx = 0
      if (.not. self%have_curr) return

      idx = find_data_index(self%data_curr, name)
   end function get_index


   subroutine clear(self)
      class(DataManager), intent(inout) :: self

      self%is_init     = .false.
      self%have_curr   = .false.
      self%have_next   = .false.
      self%dt_main     = 0_lk
      self%switch_time = huge(1_lk)
      self%y_active    = 0

      if (allocated(self%data_curr)) deallocate(self%data_curr)
      if (allocated(self%data_next)) deallocate(self%data_next)
   end subroutine clear


   subroutine set_error_mode(self, on)
      class(DataManager), intent(inout) :: self
      logical,            intent(in)    :: on

      self%stop_on_error = on
   end subroutine set_error_mode


   integer function get_sim_calendar(self) result(cal_kind)
      class(DataManager), intent(in) :: self

      if (.not. self%is_init) then
         cal_kind = cal_unknown
      else
         cal_kind = self%state%sim_cal%kind
      end if
   end function get_sim_calendar


   pure integer(lk) function compute_switch_time(cal, start_datetime, next_year) result(t_switch)
      type(CFCalendar), intent(in) :: cal
      type(DateTime),   intent(in) :: start_datetime
      integer,          intent(in) :: next_year

      type(DateTime) :: dtb

      dtb%year   = next_year
      dtb%month  = 1
      dtb%day    = 1
      dtb%hour   = 0
      dtb%minute = 0
      dtb%second = 0

      t_switch = seconds_between_datetimes(cal, dtb, start_datetime)
   end function compute_switch_time


   subroutine load_full_data(state, full, ok, errmsg)
      type(DataLoaderState), intent(in)    :: state
      type(InputData),       intent(inout) :: full
      logical,               intent(out)   :: ok
      character(*),          intent(out)   :: errmsg

      logical :: lok
      character(len=512) :: lmsg

      ok = .false.
      errmsg = ''

      if (allocated(full%vars)) deallocate(full%vars)

      call load_input_data(state, 1, full, lok, lmsg)

      if (.not. lok) then
         errmsg = 'load_full_data: full-period input data load failed: '//trim(lmsg)
         return
      end if

      ok = .true.
   end subroutine load_full_data


   subroutine append_input_data(full, part)
      type(InputData), intent(inout) :: full
      type(InputData), intent(in)    :: part

      integer :: i

      if (.not. allocated(part%vars)) return

      if (.not. allocated(full%vars)) then
         allocate(full%vars(size(part%vars)))

         do i = 1, size(part%vars)
            call copy_var_series(part%vars(i), full%vars(i))
         end do

         return
      end if

      if (size(full%vars) /= size(part%vars)) then
         error stop 'append_input_data: inconsistent number of variables'
      end if

      do i = 1, size(full%vars)
         call append_var_series(full%vars(i), part%vars(i))
      end do
   end subroutine append_input_data


   subroutine append_var_series(full, part)
      type(DataVarSeries), intent(inout) :: full
      type(DataVarSeries), intent(in)    :: part

      integer :: n0, n1, nadd, istart

      if (part%is_const) then
         if (.not. allocated(full%name)) full%name = part%name
         if (.not. allocated(full%units) .and. allocated(part%units)) full%units = part%units

         full%is_const    = .true.
         full%const_value = part%const_value
         full%n           = 0
         full%idx         = 1
         full%t_next      = huge(1_lk)

         if (allocated(full%t_axis)) deallocate(full%t_axis)
         if (allocated(full%t_edge)) deallocate(full%t_edge)
         if (allocated(full%values)) deallocate(full%values)

         return
      end if

      if (full%is_const) return

      if (.not. allocated(full%t_axis)) then
         call copy_var_series(part, full)
         return
      end if

      if (.not. allocated(part%t_axis)) return

      n0 = size(full%t_axis)
      n1 = size(part%t_axis)

      if (n1 <= 0) return

      istart = 1

      ! Avoid duplicated records if adjacent yearly slices overlap.
      if (n0 > 0 .and. part%t_axis(1) <= full%t_axis(n0)) then
         do while (istart <= n1 .and. part%t_axis(istart) <= full%t_axis(n0))
            istart = istart + 1
         end do
      end if

      if (istart > n1) return

      nadd = n1 - istart + 1

      full%t_axis = [full%t_axis, part%t_axis(istart:n1)]
      full%values = [full%values, part%values(istart:n1)]
      full%n      = n0 + nadd
   end subroutine append_var_series


   subroutine copy_var_series(src, dst)
      type(DataVarSeries), intent(in)    :: src
      type(DataVarSeries), intent(inout) :: dst

      if (allocated(src%name))  dst%name  = src%name
      if (allocated(src%units)) dst%units = src%units

      dst%is_const    = src%is_const
      dst%const_value = src%const_value
      dst%idx         = 1
      dst%n           = src%n
      dst%t_next      = src%t_next
      dst%cal         = src%cal
      dst%u           = src%u
      dst%sim_offset  = src%sim_offset
      dst%time_mode   = src%time_mode
      dst%repeat_year = src%repeat_year

      if (allocated(dst%t_axis)) deallocate(dst%t_axis)
      if (allocated(dst%t_edge)) deallocate(dst%t_edge)
      if (allocated(dst%values)) deallocate(dst%values)

      if (allocated(src%t_axis)) then
         allocate(dst%t_axis(size(src%t_axis)))
         dst%t_axis = src%t_axis
      end if

      if (allocated(src%values)) then
         allocate(dst%values(size(src%values)))
         dst%values = src%values
      end if

      if (allocated(src%t_edge)) then
         allocate(dst%t_edge(size(src%t_edge)))
         dst%t_edge = src%t_edge
      end if

      ! Edges are rebuilt after full concatenation.
   end subroutine copy_var_series

   subroutine copy_loaded_repeat_series(state, src, dst)
      type(DataLoaderState), intent(in)    :: state
      type(InputData),       intent(in)    :: src
      type(InputData),       intent(inout) :: dst

      integer :: i

      if (.not. allocated(src%vars)) return

      if (.not. allocated(dst%vars)) then
         allocate(dst%vars(size(src%vars)))
      else if (size(dst%vars) /= size(src%vars)) then
         deallocate(dst%vars)
         allocate(dst%vars(size(src%vars)))
      end if

      do i = 1, size(state%specs)
         dst%vars(i)%name = state%specs(i)%name

         if (state%specs(i)%input_type == DATA_INPUT_FILE .and. &
            state%specs(i)%repeat_enabled) then
            call copy_var_series(src%vars(i), dst%vars(i))
         end if
      end do
   end subroutine copy_loaded_repeat_series


   subroutine initialise_full_cursors(input)
      type(InputData), intent(inout) :: input

      integer :: i

      if (.not. allocated(input%vars)) return

      do i = 1, size(input%vars)
         if (.not. input%vars(i)%is_const) call rebuild_edges(input%vars(i))
         call init_series_cursor(input%vars(i))
      end do
   end subroutine initialise_full_cursors


   subroutine rebuild_edges(series)
      type(DataVarSeries), intent(inout) :: series

      integer :: i, n
      integer(lk) :: dt_first, dt_last

      if (series%is_const) return
      if (.not. allocated(series%t_axis)) return

      n = size(series%t_axis)
      series%n = n

      if (allocated(series%t_edge)) deallocate(series%t_edge)
      allocate(series%t_edge(n + 1))

      if (n >= 2) then
         dt_first = max(1_lk, series%t_axis(2) - series%t_axis(1))
         dt_last  = max(1_lk, series%t_axis(n) - series%t_axis(n-1))

         series%t_edge(1)   = series%t_axis(1) - dt_first / 2_lk
         series%t_edge(n+1) = series%t_axis(n) + dt_last  / 2_lk

         do i = 1, n - 1
            series%t_edge(i+1) = (series%t_axis(i) + series%t_axis(i+1)) / 2_lk
         end do
      else if (n == 1) then
         series%t_edge(1) = series%t_axis(1) - huge(1_lk) / 4_lk
         series%t_edge(2) = series%t_axis(1) + huge(1_lk) / 4_lk
      end if
   end subroutine rebuild_edges


   subroutine clear_input_data(input)
      type(InputData), intent(inout) :: input

      if (allocated(input%vars)) deallocate(input%vars)
   end subroutine clear_input_data


   subroutine stop_fatal(where, msg, stop_on_error)
      character(*), intent(in) :: where, msg
      logical,      intent(in) :: stop_on_error

      if (.not. stop_on_error) return

      write(*,*) 'FATAL DataManager['//trim(where)//']: ', trim(msg)
      stop 1
   end subroutine stop_fatal

end module data_manager