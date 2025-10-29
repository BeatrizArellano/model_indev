module sim_clocks
  !! Time-step management
  use precision_types,  only: rk, lk
  use time_types,       only: DateTime, CFCalendar
  use cf_time_utils,    only: seconds_between_datetimes
  implicit none
  private
  public :: init_clock, year_to_simyear

contains

  subroutine init_clock(cal, t_start, t_end, dt_main,                    &
                        sim_time_sec, n_steps, last_dt_length)
    !! Number of steps betweeen [t_start, t_end] (inclusive of start, exclusive of end).
    !!
    !! Inputs:
    !!   cal        : CFCalendar   (calendar to interpret DateTimes)
    !!   t_start    : DateTime     (simulation start)
    !!   t_end      : DateTime     (simulation end; must be > t_start)
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
    type(DateTime),   intent(in)  :: t_start, t_end
    integer(lk),      intent(in)  :: dt_main
    integer(lk),      intent(out) :: sim_time_sec, last_dt_length, n_steps

    if (dt_main <= 0_lk) error stop 'init_clock: dt must be > 0.'

    ! Convert to absolute seconds in the chosen calendar
    sim_time_sec = seconds_between_datetimes(cal, t_end, t_start)

    if (sim_time_sec <= 0_lk) error stop 'init_clock: end must be > start.'

    ! Calculates the ceiling of n_steps
    n_steps = (sim_time_sec + dt_main - 1_lk) / dt_main
    last_dt_length = sim_time_sec - (n_steps-1_lk)*dt_main

  end subroutine init_clock

  ! Gets the simulation year from a given year
    pure integer function year_to_simyear(y, start_year, end_year) result(simyear)
        integer, intent(in) :: y, start_year, end_year
        if (y < start_year .or. y > end_year) then
            simyear = -1
        else
            simyear = (y - start_year) + 1
        end if
    end function year_to_simyear

end module sim_clocks
