module sim_clocks
  !! Time-step management
  use precision_types,  only: rk, lk
  use time_types,       only: DateTime, CFCalendar, CFUnits
  use time_utils,       only: datetime_to_str
  use cf_time_utils,    only: seconds_between_datetimes, seconds_to_datetime
  use, intrinsic        :: iso_fortran_env, only: output_unit
  implicit none
  private
  public :: init_clock, year_to_simyear, simtime_to_datetime, print_progress

contains

  subroutine init_clock(cal, start_datetime, t_end, dt_main,                    &
                        sim_time_sec, n_steps, last_dt_length)
    !! Number of steps betweeen [start_datetime, t_end] (inclusive of start, exclusive of end).
    !!
    !! Inputs:
    !!   cal        : CFCalendar   (calendar to interpret DateTimes)
    !!   start_datetime    : DateTime     (simulation start)
    !!   t_end      : DateTime     (simulation end; must be > start_datetime)
    !!   dt         : int(lk)      (time-step in seconds; must be > 0)
    !!
    !! Outputs:    
    !!   sim_time_sec  : integer(lk)   Legth of simulation in seconds
    !!   n_steps       : integer         total number of time-steps
    !!   last_dt_length : integer(lk)  duration of the final step in seconds (≤ dt_sec; equals
    !!                                   dt_sec when span is an exact multiple)
    !! - Loop will typically be: do i = 1, n_steps; t = t0_sec + (i-1)*dt_sec
    !!
    type(CFCalendar), intent(in)  :: cal
    type(DateTime),   intent(in)  :: start_datetime, t_end
    integer(lk),      intent(in)  :: dt_main
    integer(lk),      intent(out) :: sim_time_sec, last_dt_length, n_steps

    if (dt_main <= 0_lk) error stop 'init_clock: dt must be > 0.'

    ! Convert to absolute seconds in the chosen calendar
    sim_time_sec = seconds_between_datetimes(cal, t_end, start_datetime)

    if (sim_time_sec <= 0_lk) error stop 'init_clock: end must be > start.'

    ! Calculates the ceiling of n_steps
    n_steps = (sim_time_sec + dt_main - 1_lk) / dt_main
    last_dt_length = sim_time_sec - (n_steps-1_lk)*dt_main

  end subroutine init_clock

  pure subroutine simtime_to_datetime(cal, start_datetime, elapsed_time, dt_out, doy)
      !! Convert elapsed_time (seconds since start_datetime) into a DateTime + doy
      type(CFCalendar), intent(in)  :: cal
      type(DateTime),   intent(in)  :: start_datetime
      integer(lk),      intent(in)  :: elapsed_time
      type(DateTime),   intent(out) :: dt_out
      integer,          intent(out) :: doy

      type(CFUnits) :: u
      integer :: y, m, d, h, mi, s

      ! Define units: "seconds since start_datetime"
      u%timeunit_to_seconds = 1.0_rk
      u%reference_year  = start_datetime%year
      u%reference_month = start_datetime%month
      u%reference_day   = start_datetime%day
      u%reference_hour  = start_datetime%hour
      u%reference_min   = start_datetime%minute
      u%reference_sec   = start_datetime%second
      u%has_time        = start_datetime%has_time

      call seconds_to_datetime(cal, u, real(elapsed_time, rk), &
                               y, m, d, h, mi, s, doy)

      dt_out%year   = y
      dt_out%month  = m
      dt_out%day    = d
      dt_out%hour   = h
      dt_out%minute = mi
      dt_out%second = s
      dt_out%has_time = .true.
  end subroutine simtime_to_datetime




  ! Gets the simulation year from a given year
    pure integer function year_to_simyear(y, start_year, end_year) result(simyear)
        integer, intent(in) :: y, start_year, end_year
        if (y < start_year .or. y > end_year) then
            simyear = -1
        else
            simyear = (y - start_year) + 1
        end if
    end function year_to_simyear

  subroutine print_progress(istep, n_steps, elapsed_time, sim_time_sec,  &
                            t_now, doy, t_start, t_end)
    integer(lk),      intent(in) :: istep, n_steps
    integer(lk),      intent(in) :: elapsed_time, sim_time_sec
    type(DateTime),   intent(in) :: t_now
    integer,          intent(in) :: doy
    type(DateTime),   intent(in) :: t_start, t_end

    real(rk) :: frac
    integer  :: pct, bar_width, filled, simyear, nyears
    character(len=16)  :: date_str      ! "YYYY-MM-DD"
    character(len=60)  :: bar
    character(len=160) :: line

    ! Persistent state to avoid spamming: last day-of-year printed
    integer, save :: last_doy = -1

    if (sim_time_sec <= 0_lk .or. n_steps <= 0_lk) return

    if (istep == 1_lk) then
       last_doy = -1
    end if

    ! Fraction completed and percentage
    if (istep == n_steps) then
        ! Force full completion on final step, even if dt does not divide the span exactly
        frac = 1.0_rk
    else
        frac = real(elapsed_time, rk) / real(sim_time_sec, rk)
        frac = max(0.0_rk, min(frac, 1.0_rk))
    end if
    pct  = int(frac*100.0_rk + 0.5_rk)

    ! Throttle to once per simulated day (unless final step)
    if (doy == last_doy .and. istep < n_steps) return
    last_doy = doy

    ! Date as YYYY-MM-DD (no time)
    write(date_str,'(I4.4,"-",I2.2,"-",I2.2)') t_now%year, t_now%month, t_now%day

    ! Simulation-year index (1..nyears)
    nyears  = t_end%year - t_start%year + 1
    simyear = year_to_simyear(t_now%year, t_start%year, t_end%year)

    bar_width = 50
    filled    = int(frac*bar_width)
    filled    = max(0, min(filled, bar_width))

    bar = repeat('#', filled) // repeat('.', bar_width - filled)

    write(line,'("[",A,"] ",I3,"%  ",A,"  Doy=",I3.3,"  Year: ",I0,"/",I0)') &
         bar(1:bar_width), pct, trim(date_str), doy, simyear, nyears

    write(output_unit,'(A,A)', advance='no') achar(13), trim(line)
    call flush(output_unit)

    if (istep == n_steps) then
       write(output_unit,*)
    end if
end subroutine print_progress





end module sim_clocks
