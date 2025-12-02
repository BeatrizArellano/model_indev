! src/physics/forcing_manager.F90
module forcing_manager
  use precision_types,     only: rk, lk
  use time_types,          only: DateTime, CFCalendar
  use read_config_yaml,    only: ConfigParams
  use geo_utils,           only: LocationInfo
  use load_forcing,        only: ForcingState, ForcingYearData, scan_and_init_forcing, &
                                 load_year_data, print_forcing_summary
  use cf_time_utils,       only: seconds_between_datetimes
  use sim_clocks,          only: year_to_simyear
  use str_utils,           only: inttostr
  implicit none
  private

  ! Data Snapshot for each model step
  type, public :: ForcingSnapshot
     real(rk) :: air_temp   = 0._rk
     real(rk) :: slp        = 0._rk
     real(rk) :: rel_hum    = 0._rk
     real(rk) :: short_rad  = 0._rk
     real(rk) :: long_rad   = 0._rk
     real(rk) :: wind_u10   = 0._rk   ! zonal (m s-1)
     real(rk) :: wind_v10   = 0._rk   ! meridional (m s-1)
     ! Freshwater fluxes [m s-1]
     real(rk) :: precip     = 0._rk   ! precipitation (downward freshwater flux)
     real(rk) :: evap       = 0._rk   ! evaporation (upward freshwater flux, negative = evaporation)
     real(rk) :: runoff     = 0._rk   ! river/runoff (downward freshwater flux)
     ! Co2
     real(rk) :: co2_air    = 0._rk
  end type ForcingSnapshot

  ! ---------------Forcing Manager-------------------------------
  type, public :: ForcingManager
     ! config / plan
     type(ForcingState) :: Fstate
     type(DateTime)     :: sim_start
     type(CFCalendar)   :: sim_cal

     ! runtime data
     type(ForcingYearData), allocatable :: Ydata_curr   ! Data arrays for current and next year
     type(ForcingYearData), allocatable :: Ydata_next   ! Both will only be alive together for a few seconds

     ! runtime state
     integer(lk) :: pad_secs    = 3600_lk    ! Numbre of seconds to load data before reaching the switch into a new year
     integer(lk) :: switch_time = huge(1_lk)
     integer     :: y_active    = 0
     logical     :: is_init     = .false.
     logical     :: have_curr   = .false.
     logical     :: have_next   = .false.
     logical     :: stop_on_error   = .true.   ! stop on error?

     ! Per-instance cushion (seconds) to de-synchronise preload across many columns (if running in parallel)
     integer(lk) :: preload_cushion = 0_lk
   contains
     procedure :: init                ! Scan and initialise forcing managers
     procedure :: prepare  
     procedure :: tick
     procedure :: sample
     procedure :: clear
     procedure :: set_error_mode
     procedure :: set_preload_cushion
     procedure :: get_sim_calendar
  end type ForcingManager

  public :: compute_switch_time  

contains

    ! -------Init: scan forcing and plan ------------------------------
    subroutine init(self, params, calendar_cfg, location, start_datetime, end_datetime, ok, errmsg)
        class(ForcingManager), intent(inout) :: self
        type(ConfigParams),    intent(in)    :: params
        type(CFCalendar),      intent(in)    :: calendar_cfg
        type(LocationInfo),    intent(in)    :: location
        type(DateTime),        intent(in)    :: start_datetime, end_datetime
        logical,               intent(out)   :: ok
        character(*),          intent(out)   :: errmsg

        ok = .false.; errmsg = ''
        write(*,'(A)') 'Scanning forcing data...'

        call scan_and_init_forcing(params, calendar_cfg, location, start_datetime, end_datetime, self%Fstate, ok, errmsg)
        if (.not. ok) then
            call stop_fatal('scan_and_init_forcing', errmsg, self%stop_on_error)
            return
        end if

        call print_forcing_summary(self%Fstate)

        self%sim_cal%kind = self%Fstate%sim_cal%kind
        self%sim_start    = start_datetime
        self%is_init      = .true.
        self%have_curr    = .false.
        self%have_next    = .false.
        ok = .true.
    end subroutine init

    ! ---------------------- Prepare: first switch and first year ------------------
    subroutine prepare(self, dt_main, preload_pad_sec, ok, errmsg)
        class(ForcingManager), intent(inout) :: self
        integer(lk),           intent(in)    :: dt_main        ! main time-step
        integer(lk), optional, intent(in)    :: preload_pad_sec
        logical,               intent(out)   :: ok
        character(*),          intent(out)   :: errmsg

        integer :: k
        logical :: lok
        character(len=256) :: lmsg

        ok = .false.; errmsg = ''
        if (.not. self%is_init) then
            errmsg = 'ForcingManager not initialized'
            call stop_fatal('prepare', errmsg, self%stop_on_error)
        end if

        ! Pad interval (ensure >= dt_main)
        if (present(preload_pad_sec)) then
            self%pad_secs = max(dt_main, preload_pad_sec)
        else
            self%pad_secs = max(dt_main, 3600_lk)
        end if

        ! Active year and switch time for next year
        self%y_active   = self%sim_start%year
        self%switch_time = compute_switch_time(self%sim_cal, self%sim_start, self%y_active + 1)

        ! Load current year
        k = year_to_simyear(self%y_active, self%Fstate%sim_y_start, self%Fstate%sim_y_end)
        if (k <= 0) then
            errmsg = 'start year outside simulation range'
            call stop_fatal('prepare: ', errmsg, self%stop_on_error)
            return
        end if

        if (.not. allocated(self%Ydata_curr)) allocate(self%Ydata_curr)
        call load_year_data(self%Fstate, k, self%Ydata_curr, lok, lmsg)
        if (.not. lok) then
            errmsg = '(Y='//trim(adjustl(inttostr(self%sim_start%year)))//') failed: '//trim(lmsg)
            call stop_fatal('load_year_data', errmsg, self%stop_on_error)
            return
        end if

        self%have_curr = .true.
        self%have_next = .false.
        ok = .true.
    end subroutine prepare

    ! ---------------------- Just in time preload and promote ----------------------
    subroutine tick(self, model_time, ok, errmsg)
        class(ForcingManager), intent(inout) :: self
        integer(lk),           intent(in)    :: model_time   ! Model time in seconds
        logical,     optional, intent(out)   :: ok
        character(*),optional, intent(out)   :: errmsg

        integer :: k
        logical :: lok
        character(len=256) :: lmsg
        integer(lk) :: guard_edge

        if (present(ok))    ok    = .true.
        if (present(errmsg)) errmsg = ''

        if (.not. self%is_init) return
        if (.not. self%have_curr) return

        ! Preload data (with optional cushion)
        guard_edge = max(0_lk, self%pad_secs + self%preload_cushion)

        if ((.not. self%have_next) .and. (model_time >= self%switch_time - guard_edge)) then
            k = year_to_simyear(self%y_active + 1, self%Fstate%sim_y_start, self%Fstate%sim_y_end)
            if (k > 0) then
                if (.not. allocated(self%Ydata_next)) allocate(self%Ydata_next)
                call load_year_data(self%Fstate, k, self%Ydata_next, lok, lmsg)
                if (.not. lok) then
                if (present(ok))     ok = .false.
                if (present(errmsg)) errmsg = 'Preload next year failed: '//trim(lmsg)
                call stop_fatal('tick', errmsg, self%stop_on_error)
                return
                end if
                self%have_next = .true.
            end if
        end if

        ! Promote Ydata_next to Ydata_curr at or after switch_time 
        if (model_time >= self%switch_time) then
            if (self%have_next) then
                if (allocated(self%Ydata_curr)) deallocate(self%Ydata_curr)
                call move_alloc(from=self%Ydata_next, to=self%Ydata_curr) ! Promotes and clears Ydata_next
                !self%Ydata_curr    = self%Ydata_next
                self%have_curr = .true.
                self%have_next = .false.
            else
                ! Load if not loaded yet
                k = year_to_simyear(self%y_active + 1, self%Fstate%sim_y_start, self%Fstate%sim_y_end)
                if (k > 0) then
                    if (.not. allocated(self%Ydata_curr)) allocate(self%Ydata_curr)
                    call load_year_data(self%Fstate, k, self%Ydata_curr, lok, lmsg)
                    if (.not. lok) then
                        if (present(ok))     ok = .false.
                        if (present(errmsg)) errmsg = 'Just in time year data load failed: '//trim(lmsg)
                        call stop_fatal('tick/just_in_time', errmsg, self%stop_on_error)
                        return
                    end if
                self%have_curr = .true.
                else
                    ! Nothing-  End of sim years
                end if
            end if

            self%y_active    = self%y_active + 1
            self%switch_time = compute_switch_time(self%sim_cal, self%sim_start, self%y_active + 1)
        end if
    end subroutine tick

    ! ------ Sampling of data per-step -----------------------
    subroutine sample(self, model_time, Fsnp, ok, errmsg)
        class(ForcingManager), intent(inout) :: self
        integer(lk),           intent(in)    :: model_time
        type(ForcingSnapshot), intent(out)   :: Fsnp
        logical,     optional, intent(out)   :: ok
        character(*),optional, intent(out)   :: errmsg

        if (present(ok))    ok    = .true.
        if (present(errmsg)) errmsg = ''

        if (.not. self%have_curr) then
            if (present(ok))     ok = .false.
            if (present(errmsg)) errmsg = 'Forcing data snapshot failed, data is not loaded'
            call stop_fatal('sample', errmsg, self%stop_on_error)
            return
        end if

        ! NOTE: value_at_step expects model_time in seconds since start of the simulaion
        Fsnp%air_temp  = self%Ydata_curr%air_temp%value_at_step(model_time)
        Fsnp%slp       = self%Ydata_curr%slp%value_at_step(model_time)
        Fsnp%rel_hum   = self%Ydata_curr%rel_hum%value_at_step(model_time)
        Fsnp%short_rad = self%Ydata_curr%short_rad%value_at_step(model_time)
        Fsnp%long_rad  = self%Ydata_curr%long_rad%value_at_step(model_time)
        Fsnp%wind_u10  = self%Ydata_curr%wind_u10%value_at_step(model_time)
        Fsnp%wind_v10  = self%Ydata_curr%wind_v10%value_at_step(model_time)
        Fsnp%precip     = self%Ydata_curr%precip%value_at_step(model_time)
        Fsnp%evap       = self%Ydata_curr%evap%value_at_step(model_time)
        Fsnp%runoff     = self%Ydata_curr%runoff%value_at_step(model_time)
        Fsnp%co2_air   = self%Ydata_curr%co2_air%value_at_step(model_time)
    end subroutine sample

    ! ------------------------------- Clear ----------------------------------
    subroutine clear(self)
        class(ForcingManager), intent(inout) :: self
        self%is_init   = .false.
        self%have_curr = .false.
        self%have_next = .false.
        if (allocated(self%Ydata_curr)) deallocate(self%Ydata_curr)
        if (allocated(self%Ydata_next)) deallocate(self%Ydata_next)
    end subroutine clear

    ! Switch to stop the full simulation in case of failure -
    ! Adding this because not sure about stopping when running parallel processes in the future
    subroutine set_error_mode(self, on)
        class(ForcingManager), intent(inout) :: self
        logical,               intent(in)    :: on
        self%stop_on_error = on
    end subroutine set_error_mode
  
    ! Sets a cushion of time to preload the data
    subroutine set_preload_cushion(self, cushion_seconds)
        class(ForcingManager), intent(inout) :: self
        integer(lk),           intent(in)    :: cushion_seconds
        self%preload_cushion = cushion_seconds
    end subroutine set_preload_cushion

    ! Returns calendar kind found in forcing data (int)
    integer function get_sim_calendar(self) result(cal_kind)
        class(ForcingManager), intent(in) :: self
        if (.not. self%is_init) then
            cal_kind = 0
        else
            cal_kind = self%Fstate%sim_cal%kind
        end if
    end function get_sim_calendar

    ! Computes the model time (in seconds) in which the next year is reached. 
    pure integer(lk) function compute_switch_time(cal, start_datetime, next_year) result(t_switch)
        type(CFCalendar), intent(in) :: cal
        type(DateTime),   intent(in) :: start_datetime
        integer,          intent(in) :: next_year
        type(DateTime) :: dtb
        dtb%year   = next_year; dtb%month = 1; dtb%day = 1
        dtb%hour   = 0; dtb%minute = 0; dtb%second = 0
        t_switch = seconds_between_datetimes(cal, dtb, start_datetime)
    end function compute_switch_time

    subroutine stop_fatal(where, msg, stop_on_error)
        character(*), intent(in) :: where, msg
        logical,      intent(in) :: stop_on_error

        if (.not. stop_on_error) return

        write(*,*) 'FATAL ForcingManager['//trim(where)//']: ', trim(msg)
        stop 1
    end subroutine stop_fatal

end module forcing_manager
