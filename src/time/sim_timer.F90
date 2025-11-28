module sim_timer
  !! Simple timers based on system_clock.
  use precision_types, only: rk, lk
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  private

  public :: Timer
  public :: timer_init, timer_start, timer_stop, timer_elapsed_s, timer_print
  public :: sim_timer_start, sim_timer_stop_and_print

  type :: Timer
     integer(lk) :: start_count = 0_lk
     integer(lk) :: end_count   = 0_lk
     integer(lk) :: count_rate  = 0_lk
     logical        :: running     = .false.
  end type Timer

  ! A dedicated timer for the whole simulation
  type(Timer), save :: sim_t

contains

  subroutine timer_init(t)
    type(Timer), intent(inout) :: t
    integer(lk) :: rate

    call system_clock(count_rate=rate)
    if (rate <= 0_lk) then
       t%count_rate = 0_lk
    else
       t%count_rate = rate
    end if
    t%start_count = 0_lk
    t%end_count   = 0_lk
    t%running     = .false.
  end subroutine timer_init


  subroutine timer_start(t)
    type(Timer), intent(inout) :: t
    integer(lk) :: count, rate

    ! Ensure we have a rate
    if (t%count_rate <= 0_lk) then
       call system_clock(count_rate=rate)
       if (rate > 0_lk) t%count_rate = rate
    end if

    call system_clock(count=count)
    t%start_count = count
    t%running     = .true.
  end subroutine timer_start


  subroutine timer_stop(t)
    type(Timer), intent(inout) :: t
    integer(lk) :: count

    call system_clock(count=count)
    t%end_count = count
    t%running   = .false.
  end subroutine timer_stop


  function timer_elapsed_s(t) result(seconds)
    !! Elapsed wall-clock time in seconds.
    type(Timer), intent(in) :: t
    real(rk) :: seconds
    integer(lk) :: c_start, c_end, rate

    rate = t%count_rate
    if (rate <= 0_lk) then
       seconds = -1.0_rk
       return
    end if

    if (t%running) then
       call system_clock(count=c_end)
       c_start = t%start_count
    else
       c_start = t%start_count
       c_end   = t%end_count
    end if

    if (c_end < c_start) then
       seconds = -1.0_rk
    else
       seconds = real(c_end - c_start, rk) / real(rate, rk)
    end if
  end function timer_elapsed_s


  subroutine timer_print(t, label)
    !! Print summary
    type(Timer), intent(in) :: t
    character(len=*), intent(in), optional :: label

    real(rk) :: seconds
    integer  :: h, m
    real(rk) :: s_rem
    character(len=64) :: lab

    seconds = timer_elapsed_s(t)
    if (seconds < 0.0_rk) then
       if (present(label)) then
          write(output_unit,*) trim(label)//': timer not available.'
       else
          write(output_unit,*) 'Timer: timer not available.'
       end if
       return
    end if

    h = int(seconds / 3600.0_rk)
    m = int( (seconds - 3600.0_rk*h) / 60.0_rk )
    s_rem = seconds - 3600.0_rk*h - 60.0_rk*m

    if (present(label)) then
       lab = trim(label)
    else
       lab = 'Elapsed time: '
    end if

    write(output_unit,'(A,": ",I2.2,"h:",I2.2,"m:",F6.3,"s (",F8.3," s)")') &
         trim(lab), h, m, s_rem, seconds
  end subroutine timer_print

  ! Interfaces
  subroutine sim_timer_start()
    !! Start timing the whole simulation.
    call timer_init(sim_t)
    call timer_start(sim_t)
  end subroutine sim_timer_start


  subroutine sim_timer_stop_and_print()
    !! Stop timing and print a summary.
    call timer_stop(sim_t)
    call timer_print(sim_t, 'Total simulation time ')
  end subroutine sim_timer_stop_and_print

end module sim_timer
