module output_manager
  use precision_types,   only: rk, lk
  use read_config_yaml,  only: ConfigParams
  use geo_utils,         only: LocationInfo
  use time_utils,        only: parse_interval_to_seconds, sec_per_hour, sec_per_day
  use grids,             only: VerticalGrid
  use output_writer,     only: OutputWriter
  use variable_registry, only: VarMetadata
  implicit none
  private
  public :: OutputManager

  ! ---- Configuration parameters
  type :: OutputConfig
     character(:), allocatable :: file          ! base filename (no .nc)
     logical                   :: overwrite     ! Whether to overwrite the file or not
     character(:), allocatable :: frequency     ! 'hourly'|'daily'|'custom'
     character(:), allocatable :: statistic     ! 'mean'|'instant'
     character(:), allocatable :: interval_txt  ! e.g. '12h' (only for custom)
     integer(lk)               :: interval_s = 0_lk
     integer                   :: sync_frequency_d = 0     ! Sync frequency in days
     integer(lk)               :: sync_frequency_s = 0_lk  ! Sync frequency in seconds: 0 is no periodic sync
  end type OutputConfig

  ! ---- Scheduler
  type :: OutputScheduler
     integer(lk) :: interval_s  = 0_lk   ! Target output interval length [s]
     integer(lk) :: model_time  = 0_lk   ! Start time of current window [s since start]
     integer(lk) :: acc_time_s  = 0_lk   ! Accumulated duration in current window [s]
  contains
     procedure :: init          => sched_init
     procedure :: advance       => sched_advance      ! returns true when window closes
     procedure :: emit_time_s   => sched_emit_time_s  ! timestamp to use (center or end)
     procedure :: reset_window  => sched_reset_window
     procedure :: window_len_s  => sched_window_len_s
  end type OutputScheduler

  ! Stores data for all variables
  type :: OutputData
     character(:), allocatable :: statistic          ! mean or instant
     integer(lk) :: steps    = 0_lk
     integer(lk) :: N        = 0_lk
     integer(lk) :: Ni       = 0_lk
     integer     :: n_centre = 0
     integer     :: n_iface  = 0
     integer     :: n_scalar = 0

     ! Sums (only used if statistic is 'mean')
     real(rk), allocatable :: centre_sum(:,:)   ! (N,  n_centre)
     real(rk), allocatable :: iface_sum(:,:)    ! (Ni, n_iface)
     real(rk), allocatable :: scalar_sum(:)     ! (n_scalar)

     ! Last values (used if statistic is 'instant')
     real(rk), allocatable :: centre_last(:,:)  ! (N,  n_centre)
     real(rk), allocatable :: iface_last(:,:)   ! (Ni, n_iface)
     real(rk), allocatable :: scalar_last(:)    ! (n_scalar)
  contains
     procedure :: init             => outdata_init
     procedure :: update_from_vars => outdata_update_from_vars
     procedure :: finalize         => outdata_finalize
     procedure :: reset            => outdata_reset
  end type OutputData

  ! --- Output Orchestrator
  type :: OutputManager
     logical :: is_init = .false.

     type(OutputConfig)    :: cfg
     type(OutputScheduler) :: sched
     type(OutputData)      :: acc
     type(OutputWriter)    :: writer

     integer(lk)          :: N  = 0_lk
     integer(lk)          :: Ni = 0_lk
     logical              :: is_mean = .true.
     integer(lk)          :: next_sync_time_s = -1_lk

     ! Variable metadata 
     type(VarMetadata), allocatable :: vars(:)
     integer :: n_centre = 0                   ! Number of variables defined at layers' centres
     integer :: n_iface  = 0                   ! Number of variables defined at layers' interfaces
     integer :: n_scalar = 0
     ! Persistent record arrats (surface to bottom orientation)
     real(rk), allocatable :: centre_rec(:,:)   ! (N,  n_centre)
     real(rk), allocatable :: iface_rec(:,:)    ! (Ni, n_iface)
     real(rk), allocatable :: scalar_rec(:)     ! (n_scalar)
  contains
     procedure :: init   => om_init      ! reads config, sets scheduler+acc, opens NetCDF
     procedure :: step   => om_step      ! update accumulators; write when window closes
     procedure :: close  => om_close     ! close NetCDF
  end type OutputManager

contains
  !================ CONFIGuration parameters ===
  subroutine read_output_config(params, cfg)
    type(ConfigParams), intent(in)  :: params
    type(OutputConfig), intent(out) :: cfg

    character(:), allocatable :: fname, frequency, statistic, interval
    logical                   :: overwrite
    character(len=6) :: freq_choices(3) = ['hourly','daily ','custom']
    character(len=7) :: stat_choices(2) = ['mean   ','instant']

    integer(lk) :: interval_s

    fname = params%get_param_str('output.file', default='output')             ! Filename
    overwrite = params%get_param_logical('output.overwrite', default=.false.) 

    frequency  = params%get_param_str('output.frequency', default='daily', &
                                      choices=freq_choices, trim_value=.true., match_case=.false.)

    statistic = params%get_param_str('output.statistic', default='mean', &
                                     choices=stat_choices, trim_value=.true., match_case=.false.)

    select case (frequency)
      case ('hourly')
         interval_s = sec_per_hour
      case ('daily')
         interval_s = sec_per_day
      case ('custom')
         interval = params%get_param_str('output.interval', default='12h', trim_value=.true.)
         call parse_interval_to_seconds(interval, interval_s)   ! e.g., "12h" -> 43200
         if (interval_s <= 0_lk) then
            error stop 'output.interval is invalid; use forms like 900s, 15m, 12h, 2d'
         end if
         cfg%interval_txt = trim(interval)
    end select

    cfg%sync_frequency_d = params%get_param_int('output.sync_frequency',default=30, min=0)
    cfg%sync_frequency_s = int(cfg%sync_frequency_d, lk) * sec_per_day

    cfg%file       = fname
    cfg%overwrite  = overwrite
    cfg%frequency  = trim(frequency)
    cfg%statistic  = trim(statistic)
    cfg%interval_s = interval_s
  end subroutine read_output_config

  !================ SCHEDULER ================
  ! Initialise the scheduler with model timestep, output interval, and starting time.
   subroutine sched_init(this, dt_s, interval_s, model_time)
      class(OutputScheduler), intent(inout) :: this
      integer(lk),           intent(in)     :: dt_s        ! <- ignored, kept for backwards compat
      integer(lk),           intent(in)     :: interval_s
      integer(lk),           intent(in), optional :: model_time

      if (interval_s <= 0_lk) error stop 'Scheduler:init interval_s must be > 0.'

      this%interval_s = interval_s
      this%model_time = merge(model_time, 0_lk, present(model_time))
      this%acc_time_s = 0_lk
   end subroutine sched_init


   ! Advance the scheduler one model time-step. Returns true when an output interval is complete.
   logical function sched_advance(this, dt_s) result(emit)
      class(OutputScheduler), intent(inout) :: this
      integer(lk),           intent(in)     :: dt_s

      if (dt_s <= 0_lk) error stop 'Scheduler:advance dt_s must be > 0.'

      this%acc_time_s = this%acc_time_s + dt_s
      emit = (this%acc_time_s >= this%interval_s)
   end function sched_advance


   ! Return the time (centre or end) of the current completed window.
   integer(lk) function sched_emit_time_s(this, is_mean) result(t_s)
      class(OutputScheduler), intent(in) :: this
      logical,                intent(in) :: is_mean
      integer(lk) :: eff_len

      if (this%acc_time_s <= 0_lk) error stop 'Scheduler: emit_time before window close.'

      eff_len = min(this%acc_time_s, this%interval_s)

      if (is_mean) then
         t_s = this%model_time + eff_len/2_lk
      else
         t_s = this%model_time + eff_len
      end if
   end function sched_emit_time_s


    ! Advance the window by one interval and reset the step counter.
   subroutine sched_reset_window(this)
      class(OutputScheduler), intent(inout) :: this
      this%model_time = this%model_time + this%acc_time_s
      this%acc_time_s = 0_lk
   end subroutine sched_reset_window


  ! Return the length of the output window in seconds.
   integer(lk) function sched_window_len_s(this) result(wlen)
      class(OutputScheduler), intent(in) :: this
      wlen = min(this%acc_time_s, this%interval_s)
   end function sched_window_len_s


  !================ Data storage for each interval ================
  ! Initialise data arrays for all varables
  subroutine outdata_init(this, N, Ni, statistic, n_centre, n_iface, n_scalar)
    class(OutputData), intent(inout) :: this
    integer(lk),       intent(in)    :: N, Ni
    character(*),      intent(in)    :: statistic
    integer,           intent(in)    :: n_centre, n_iface, n_scalar

    this%statistic = trim(statistic)
    this%steps     = 0_lk
    this%N         = N
    this%Ni        = Ni
    this%n_centre  = n_centre
    this%n_iface   = n_iface
    this%n_scalar  = n_scalar

    ! Deallocate any previous content
    if (allocated(this%centre_sum))  deallocate(this%centre_sum)
    if (allocated(this%iface_sum))   deallocate(this%iface_sum)
    if (allocated(this%scalar_sum))  deallocate(this%scalar_sum)
    if (allocated(this%centre_last)) deallocate(this%centre_last)
    if (allocated(this%iface_last))  deallocate(this%iface_last)
    if (allocated(this%scalar_last)) deallocate(this%scalar_last)

   if (this%statistic == 'mean') then
      ! Allocate only sums
      if (n_centre > 0) then
         allocate(this%centre_sum(N, n_centre))
         this%centre_sum = 0.0_rk
      end if
      if (n_iface > 0) then
         allocate(this%iface_sum(Ni, n_iface))
         this%iface_sum = 0.0_rk
      end if
      if (n_scalar > 0) then
         allocate(this%scalar_sum(n_scalar))
         this%scalar_sum = 0.0_rk
      end if
   else
      ! 'instant': allocate only last
      if (n_centre > 0) then
         allocate(this%centre_last(N, n_centre))
         this%centre_last = 0.0_rk
      end if
      if (n_iface > 0) then
         allocate(this%iface_last(Ni, n_iface))
         this%iface_last = 0.0_rk
      end if
      if (n_scalar > 0) then
         allocate(this%scalar_last(n_scalar))
         this%scalar_last = 0.0_rk
      end if
   end if
  end subroutine outdata_init

   ! Add a new model time step to the data arrays or update last values.
  subroutine outdata_update_from_vars(this, vars)
    class(OutputData),   intent(inout) :: this
    type(VarMetadata),   intent(in)    :: vars(:)

    integer :: j, ic, ii, is

    ic = 0; ii = 0; is = 0

    do j = 1, size(vars)
       if (.not. vars(j)%output) cycle

       select case (vars(j)%n_space_dims)
       case (1)
          select case (trim(vars(j)%vert_coord))
          case ('centre')
             ic = ic + 1
             if (.not. associated(vars(j)%data_1d)) error stop 'OutputData:update: centre var has no data_1d'
             if (allocated(this%centre_sum)) then
                this%centre_sum(:, ic) = this%centre_sum(:, ic) + vars(j)%data_1d(:)
             else
                this%centre_last(:, ic) = vars(j)%data_1d(:)
             end if

          case ('interface')
             ii = ii + 1
             if (.not. associated(vars(j)%data_1d)) error stop 'OutputData:update: iface var has no data_1d'
             if (allocated(this%iface_sum)) then
                this%iface_sum(:, ii) = this%iface_sum(:, ii) + vars(j)%data_1d(:)
             else
                this%iface_last(:, ii) = vars(j)%data_1d(:)
             end if

          case default
             error stop 'OutputData:update: unsupported vert_coord for profile.'
          end select

       case (0)
          is = is + 1
          if (.not. associated(vars(j)%data_0d)) error stop 'OutputData:update: scalar var has no data_0d'
          if (allocated(this%scalar_sum)) then
             this%scalar_sum(is) = this%scalar_sum(is) + vars(j)%data_0d
          else
             this%scalar_last(is) = vars(j)%data_0d
          end if

       case default
          error stop 'OutputData:update: unsupported n_space_dims.'
       end select
    end do

    this%steps = this%steps + 1_lk
  end subroutine outdata_update_from_vars

  ! Produce the final interval values (mean or last) for all variables
  subroutine outdata_finalize(this, centre_out, iface_out, scalar_out)
    class(OutputData), intent(in)  :: this
    real(rk),          intent(out), optional :: centre_out(:,:)
    real(rk),          intent(out), optional :: iface_out(:,:)
    real(rk),          intent(out), optional :: scalar_out(:)

    real(rk) :: denom

    if (this%steps <= 0_lk) error stop 'OutputData: finalize empty window'

    if (this%statistic == 'mean') then
       denom = real(max(1_lk, this%steps), rk)
    else
       denom = 1.0_rk
    end if

    if (present(centre_out) .and. this%n_centre > 0) then
       if (allocated(this%centre_sum)) then
          centre_out = this%centre_sum / denom
       else
          centre_out = this%centre_last
       end if
    end if

    if (present(iface_out) .and. this%n_iface > 0) then
       if (allocated(this%iface_sum)) then
          iface_out = this%iface_sum / denom
       else
          iface_out = this%iface_last
       end if
    end if

    if (present(scalar_out) .and. this%n_scalar > 0) then
       if (allocated(this%scalar_sum)) then
          scalar_out = this%scalar_sum / denom
       else
          scalar_out = this%scalar_last
       end if
    end if
  end subroutine outdata_finalize
  ! Reset accumulated sums and counters at the start of a new output interval.
  subroutine outdata_reset(this)
    class(OutputData), intent(inout) :: this
    this%steps = 0_lk
    if (allocated(this%centre_sum)) this%centre_sum = 0.0_rk
    if (allocated(this%iface_sum))  this%iface_sum  = 0.0_rk
    if (allocated(this%scalar_sum)) this%scalar_sum = 0.0_rk
    if (allocated(this%centre_last)) this%centre_last = 0.0_rk
    if (allocated(this%iface_last))  this%iface_last  = 0.0_rk
    if (allocated(this%scalar_last)) this%scalar_last = 0.0_rk
  end subroutine outdata_reset

  !================ OUTPUT MANAGER ================
  ! Read output config, set up scheduler and data storage, and open the NetCDF file.
  subroutine om_init(this, params, grid_phys, dt_s, time_units, calendar_name, vars, loc)
    class(OutputManager), intent(inout) :: this
    type(ConfigParams),   intent(in)    :: params
    type(VerticalGrid),   intent(in)    :: grid_phys
    integer(lk),          intent(in)    :: dt_s
    character(*),         intent(in)    :: time_units      ! "seconds since [datetime]"
    character(*),         intent(in)    :: calendar_name
    type(VarMetadata),    intent(in)    :: vars(:)
    type(LocationInfo),   intent(in)    :: loc

    character(:), allocatable :: path
    character(:), allocatable :: title_str
    logical :: exists
    integer :: ios
    integer :: j

    ! Read config
    call read_output_config(params, this%cfg)

    ! Build full output path (base name + .nc)
    path = trim(this%cfg%file)//'.nc'

    inquire(file=trim(path), exist=exists)
    if (exists) then
       if (.not. this%cfg%overwrite) then
          write(*,*) 'ERROR: output file ', trim(path), ' already exists and overwrite = no.'
          stop 1
       else
          write(*,*) 'Overwriting existing output file: ', trim(path)
          call delete_file_if_exists(path, ios)
          if (ios /= 0) then
             write(*,*) 'ERROR: Failed to delete old output file: ', trim(path), ' iostat=', ios
             stop 1
          end if
       end if
    else
       write(*,*) 'Creating new output file: ', trim(path)
    end if

    this%N      = grid_phys%nz
    this%Ni     = grid_phys%nz + 1
    this%is_mean = (this%cfg%statistic == 'mean')

    ! Store metadata (shallow copy: pointers stay pointing to live arrays)
    if (allocated(this%vars)) deallocate(this%vars)
    allocate(this%vars(size(vars)))
    this%vars = vars

    ! Count how many of each type we have (centre/interface/scalar)
    this%n_centre = 0
    this%n_iface  = 0
    this%n_scalar = 0

    do j = 1, size(this%vars)
       if (.not. this%vars(j)%output) cycle

       select case (this%vars(j)%n_space_dims)
       case (1)
          select case (trim(this%vars(j)%vert_coord))
          case ('centre')
             this%n_centre = this%n_centre + 1
          case ('interface')
             this%n_iface  = this%n_iface  + 1
          case default
             error stop 'OutputManager:init: unsupported vert_coord for profile.'
          end select
       case (0)
          this%n_scalar = this%n_scalar + 1
       case default
          error stop 'OutputManager:init: unsupported n_space_dims.'
       end select
    end do

    ! Init scheduler and arrays for sums (for mean statistic)
    call this%sched%init(dt_s, this%cfg%interval_s, model_time=0_lk)
    call this%acc%init(N       = this%N,         &
                       Ni      = this%Ni,        &
                       statistic = this%cfg%statistic, &
                       n_centre = this%n_centre, &
                       n_iface  = this%n_iface,  &
                       n_scalar = this%n_scalar)

    ! Allocate persistent record arrays
    if (allocated(this%centre_rec)) deallocate(this%centre_rec)
    if (allocated(this%iface_rec))  deallocate(this%iface_rec)
    if (allocated(this%scalar_rec)) deallocate(this%scalar_rec)

    if (this%n_centre > 0) allocate(this%centre_rec(this%N,  this%n_centre))
    if (this%n_iface  > 0) allocate(this%iface_rec (this%Ni, this%n_iface))
    if (this%n_scalar > 0) allocate(this%scalar_rec(this%n_scalar))

    ! Periodic sync if requested
    if (this%cfg%sync_frequency_s > 0_lk) then
       this%next_sync_time_s = this%cfg%sync_frequency_s
    else
       this%next_sync_time_s = -1_lk
    end if

    call build_output_title(loc, title_str)

    ! Open NetCDF file and define structure using the same vars
    call this%writer%open_file(path              = path,              &
                               grid              = grid_phys,         &
                               interval_statistic= this%cfg%statistic,&
                               time_units        = time_units,        &
                               calendar_name     = calendar_name,     &
                               vars              = this%vars,         &
                               title             = title_str)

    this%is_init = .true.
  end subroutine om_init

  ! Update data storage arrays each timestep and write a record when an output interval closes.
   subroutine om_step(this, dt_s)
      class(OutputManager), intent(inout) :: this
      integer(lk),          intent(in)    :: dt_s

      logical :: emit
      integer(lk) :: t_emit, wlen

      if (.not. this%is_init) error stop 'OutputManager: step called before init'
      if (dt_s <= 0_lk) error stop 'OutputManager: step dt_s must be > 0.'

      ! Update accumulators from current instantaneous state
      call this%acc%update_from_vars(this%vars)

      ! Advance scheduler with actual dt for this step
      emit = this%sched%advance(dt_s)

      if (emit) then
         
         ! Finalize statistics for current window
         call this%acc%finalize(centre_out = this%centre_rec, &
                                 iface_out  = this%iface_rec,  &
                                 scalar_out = this%scalar_rec)

         ! Time info for this record
         t_emit = this%sched%emit_time_s(this%is_mean)
         wlen   = this%sched%window_len_s()

         ! Append to file
         call this%writer%append_record(elapsed_time_s = t_emit,        &
                                       window_len_s   = wlen,          &
                                       is_mean        = this%is_mean,  &
                                       centre_data    = this%centre_rec, &
                                       iface_data     = this%iface_rec,  &
                                       scalar_data    = this%scalar_rec)

         ! Optional periodic sync to disk
         if (this%next_sync_time_s > 0_lk .and. t_emit >= this%next_sync_time_s) then
            call this%writer%sync_file()
            this%next_sync_time_s = this%next_sync_time_s + this%cfg%sync_frequency_s
         end if

         ! Reset for next window
         call this%acc%reset()
         call this%sched%reset_window()
      end if
   end subroutine om_step


  ! Finalise and close the NetCDF output file (optionally syncing to disk).
  subroutine om_close(this, sync_now)
    class(OutputManager), intent(inout) :: this
    logical, intent(in), optional :: sync_now
    if (.not. this%is_init) return
    call this%writer%close_file(do_sync=merge(sync_now, .false., present(sync_now)))
    this%is_init = .false.
  end subroutine om_close

  ! Delete an existing file if it is present on disk.
  subroutine delete_file_if_exists(path, iostat_out)
    character(*), intent(in)  :: path
    integer,      intent(out) :: iostat_out

    logical :: exists
    integer :: u, ios

    iostat_out = 0
    inquire(file=trim(path), exist=exists)
    if (.not. exists) return

    open(newunit=u, file=trim(path), status='old', iostat=ios)
    if (ios /= 0) then
       iostat_out = ios
       return
    end if

    close(u, status='delete', iostat=ios)
    if (ios /= 0) then
       iostat_out = ios
    end if
  end subroutine delete_file_if_exists

    ! Build a title from LocationInfo
   subroutine build_output_title(loc, title)
      type(LocationInfo), intent(in)  :: loc
      character(:), allocatable, intent(out) :: title

      character(len=160) :: buf
      logical :: has_name

      has_name = allocated(loc%name) .and. len_trim(loc%name) > 0

      if (has_name) then
         write(buf,'(A," (lat=",F7.3,", lon=",F8.3,", depth=",F7.1," m)")') &
               trim(loc%name), loc%lat, loc%lon, loc%depth
      else
         write(buf,'("1D simulation (lat=",F7.3,", lon=",F8.3,", depth=",F7.1," m)")') &
               loc%lat, loc%lon, loc%depth
      end if

      title = trim(buf)
   end subroutine build_output_title


end module output_manager
