module output_manager
  use precision_types,  only: rk, lk
  use read_config_yaml, only: ConfigParams
  use time_utils,       only: parse_interval_to_seconds, sec_per_hour, sec_per_day
  use grids,            only: VerticalGrid
  use output_writer,    only: OutputWriter
  implicit none
  private
  public :: OutputManager


  ! ---------- Parsed config ----------
  type :: OutputConfig
     character(:), allocatable :: file          ! base filename (no .nc)
     character(:), allocatable :: frequency     ! 'hourly'|'daily'|'custom'
     character(:), allocatable :: statistic     ! 'mean'|'instant'
     character(:), allocatable :: interval_txt  ! e.g. '12h' (only for custom)
     integer(lk)               :: interval_s = 0_lk
  end type OutputConfig

  ! ---------- Scheduler (time windows only) ----------
  type :: OutputScheduler
     integer(lk) :: dt_s = 0_lk
     integer(lk) :: interval_s = 0_lk
     integer(lk) :: interval_steps = 0_lk
     integer(lk) :: steps_in_window = 0_lk
     integer(lk) :: t0_s = 0_lk
  contains
     procedure :: init          => sched_init
     procedure :: advance       => sched_advance      ! returns .true. when window closes
     procedure :: emit_time_s   => sched_emit_time_s  ! timestamp to use (center or end)
     procedure :: reset_window  => sched_reset_window
     procedure :: window_len_s  => sched_window_len_s
  end type OutputScheduler

  ! ---------- Accumulator (two diagnostics for now) ----------
  type :: OutputData
     character(:), allocatable :: statistic          ! 'mean'|'instant'
     integer(lk)               :: steps = 0_lk
     real(rk), allocatable     :: temp_sum(:), kz_sum(:)
     real(rk), allocatable     :: temp_last(:), kz_last(:)
  contains
     procedure :: init      => outdata_init
     procedure :: update    => outdata_update
     procedure :: finalize  => outdata_finalize
     procedure :: reset     => outdata_reset
  end type OutputData

  ! ---------- Orchestrator ----------
  type :: OutputManager
     logical :: is_init = .false.
     type(OutputConfig)   :: cfg
     type(OutputScheduler):: sched
     type(OutputData)     :: acc
     type(OutputWriter)   :: writer
     integer(lk)          :: N  = 0_lk
     integer(lk)          :: Ni = 0_lk
     logical              :: is_mean = .true.
  contains
     procedure :: init   => om_init      ! reads config, sets scheduler+acc, opens NetCDF
     procedure :: step   => om_step      ! update accumulators; write when window closes
     procedure :: close  => om_close     ! close NetCDF
  end type OutputManager

contains

    subroutine read_output_config(params, cfg)
        type(ConfigParams), intent(in) :: params
        type(OutputConfig), intent(out) :: cfg

        character(:), allocatable :: fname, frequency, statistic, interval
        character(len=6) :: freq_choices(3) = ['hourly','daily ','custom']
        character(len=7) :: stat_choices(2) = ['mean   ','instant']
     
        integer(lk) :: interval_s

        fname = params%get_param_str('output.file', default='output')
        frequency  = params%get_param_str('output.frequency', default='daily', &
                                           choices=freq_choices, trim_value=.true., match_case=.false.)

        statistic = params%get_param_str('output.statistic', default='mean', &
                                          choices=stat_choices, trim_value=.true., match_case=.false.)

        select case(frequency)
            case('hourly') 
                interval_s = sec_per_hour
            case('daily')
                interval_s = sec_per_day
            case('custom')
                interval = params%get_param_str('output.interval', default='12h', trim_value=.true.)
                call parse_interval_to_seconds(interval, interval_s)   ! e.g., "12h" -> 43200
                if (interval_s <= 0_lk) then
                    error stop 'output.interval is invalid; use forms like 900s, 15m, 12h, 2d'
                end if
                cfg%interval_txt = trim(interval)
        end select
        cfg%file         = fname
        cfg%frequency    = trim(frequency)
        cfg%statistic    = trim(statistic)       
        cfg%interval_s   = interval_s
    end subroutine read_output_config

    !================ SCHEDULER ================
  subroutine sched_init(this, dt_s, interval_s, t0_s)
    class(OutputScheduler), intent(inout) :: this
    integer(lk), intent(in) :: dt_s, interval_s
    integer(lk), intent(in), optional :: t0_s
    if (interval_s<=0_lk .or. dt_s<=0_lk) error stop 'Scheduler:init bad dt/interval'
    if (mod(interval_s, dt_s) /= 0_lk)     error stop 'Scheduler:init require interval % dt == 0'
    this%dt_s = dt_s
    this%interval_s = interval_s
    this%interval_steps = interval_s / dt_s
    this%steps_in_window = 0_lk
    this%t0_s = merge(t0_s, 0_lk, present(t0_s))
  end subroutine sched_init

  logical function sched_advance(this, dt_s) result(emit)
    class(OutputScheduler), intent(inout) :: this
    integer(lk), intent(in) :: dt_s
    if (dt_s /= this%dt_s) error stop 'Scheduler:advance dt mismatch'
    this%steps_in_window = this%steps_in_window + 1_lk
    emit = (this%steps_in_window == this%interval_steps)
  end function sched_advance

  integer(lk) function sched_emit_time_s(this, is_mean) result(t_s)
    class(OutputScheduler), intent(in) :: this
    logical, intent(in) :: is_mean
    if (this%steps_in_window /= this%interval_steps) error stop 'Scheduler: emit_time before window close'
    if (is_mean) then
      t_s = this%t0_s + this%interval_s/2_lk
    else
      t_s = this%t0_s + this%interval_s
    end if
  end function sched_emit_time_s

  subroutine sched_reset_window(this)
    class(OutputScheduler), intent(inout) :: this
    this%t0_s = this%t0_s + this%interval_s
    this%steps_in_window = 0_lk
  end subroutine sched_reset_window

  integer(lk) function sched_window_len_s(this) result(wlen)
    class(OutputScheduler), intent(in) :: this
    wlen = this%interval_s
  end function sched_window_len_s

  !================ ACCUMULATOR ================
  subroutine outdata_init(this, N, Ni, statistic)
    class(OutputData), intent(inout) :: this
    integer(lk), intent(in) :: N, Ni
    character(*), intent(in):: statistic
    this%statistic = trim(statistic)
    this%steps = 0_lk
    allocate(this%temp_last(N), this%kz_last(Ni))
    this%temp_last = 0.0_rk; this%kz_last = 0.0_rk
    if (this%statistic == 'mean') then
      allocate(this%temp_sum(N), this%kz_sum(Ni))
      this%temp_sum = 0.0_rk;   this%kz_sum = 0.0_rk
    else
      if (allocated(this%temp_sum)) deallocate(this%temp_sum)
      if (allocated(this%kz_sum))   deallocate(this%kz_sum)
    end if
  end subroutine outdata_init

  subroutine outdata_update(this, temp, kz)
    class(OutputData), intent(inout) :: this
    real(rk), intent(in) :: temp(:), kz(:)
    if (allocated(this%temp_sum)) then
      this%temp_sum = this%temp_sum + temp
      this%kz_sum   = this%kz_sum   + kz
    else
      this%temp_last = temp
      this%kz_last   = kz
    end if
    this%steps = this%steps + 1_lk
  end subroutine outdata_update

  subroutine outdata_finalize(this, temp_out, kz_out)
    class(OutputData), intent(in)  :: this
    real(rk), intent(out) :: temp_out(:), kz_out(:)
    if (this%steps <= 0_lk) error stop 'OutputData: finalize empty window'
    if (allocated(this%temp_sum)) then
      temp_out = this%temp_sum / real(max(1_lk,this%steps), rk)
      kz_out   = this%kz_sum   / real(max(1_lk,this%steps), rk)
    else
      temp_out = this%temp_last
      kz_out   = this%kz_last
    end if
  end subroutine outdata_finalize

  subroutine outdata_reset(this)
    class(OutputData), intent(inout) :: this
    this%steps = 0_lk
    if (allocated(this%temp_sum)) then
      this%temp_sum = 0.0_rk
      this%kz_sum   = 0.0_rk
    end if
  end subroutine outdata_reset

  !================ OUTPUT MANAGER ================
  subroutine om_init(this, params, grid_phys, dt_s, time_units, calendar_name, title)
    class(OutputManager), intent(inout) :: this
    type(ConfigParams),   intent(in)    :: params
    type(VerticalGrid),   intent(in)    :: grid_phys
    integer(lk),          intent(in)    :: dt_s
    character(*),         intent(in)    :: time_units      ! "seconds since ..."
    character(*),         intent(in)    :: calendar_name
    character(*),         intent(in), optional :: title

    call read_output_config(params, this%cfg)
    this%N  = grid_phys%nz
    this%Ni = grid_phys%nz + 1
    this%is_mean = (this%cfg%statistic == 'mean')

    call this%sched%init(dt_s, this%cfg%interval_s, t0_s=0_lk)
    call this%acc%init(N=this%N, Ni=this%Ni, statistic=this%cfg%statistic)

    call this%writer%open_file(path=trim(this%cfg%file)//'.nc', grid_phys=grid_phys, &
                               interval_statistic=this%cfg%statistic, time_units=time_units, &
                               calendar_name=calendar_name, title=merge(title,'',present(title)))
    this%is_init = .true.
  end subroutine om_init

  subroutine om_step(this, temp_phys, kz_phys)
    class(OutputManager), intent(inout) :: this
    real(rk), intent(in) :: temp_phys(:), kz_phys(:)  ! bottom→surface (N), (N+1)
    real(rk), allocatable :: temp_rec(:), kz_rec(:)
    logical :: emit
    if (.not. this%is_init) error stop 'OutputManager: step called before init'

    call this%acc%update(temp_phys, kz_phys)
    emit = this%sched%advance(this%sched%dt_s)

    if (emit) then
       allocate(temp_rec(this%N), kz_rec(this%Ni))
       call this%acc%finalize(temp_rec, kz_rec)

       call this%writer%append_record( t_window_start_s = this%sched%t0_s,          &
                                       window_len_s     = this%sched%window_len_s(),&
                                       temp_phys        = temp_rec,                 &
                                       kz_phys          = kz_rec,                   &
                                       is_mean          = this%is_mean )

       call this%acc%reset()
       call this%sched%reset_window()
    end if
  end subroutine om_step

  subroutine om_close(this, sync_now)
    class(OutputManager), intent(inout) :: this
    logical, intent(in), optional :: sync_now
    if (.not. this%is_init) return
    call this%writer%close_file(do_sync=merge(sync_now,.false.,present(sync_now)))
    this%is_init = .false.
  end subroutine om_close



end module output_manager
