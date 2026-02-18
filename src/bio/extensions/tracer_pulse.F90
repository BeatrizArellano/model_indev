module tracer_pulse
   use cf_time_utils,    only: seconds_between_datetimes
   use fabm,             only: type_fabm_model
   use find_utils,       only: find_nearest_index
   use grids,            only: VerticalGrid
   use precision_types,  only: rk, lk
   use read_config_yaml, only: ConfigParams
   use str_utils,        only: to_lower, realtostr
   use time_utils,       only: datetime_from_string, is_datetime_before, is_datetime_equal, &
                               parse_interval_to_seconds
   use time_types,       only: DateTime, CFCalendar

   implicit none
   private
   public :: tracer_pulse_init, tracer_pulse_apply, tracer_pulse_clear
   public :: TracerPulse

    type :: TracerPulse
        ! Tracer
        integer :: tracer_idx = -1
        character(len=:), allocatable :: tracer_name
        ! Time
        type(DateTime) :: pulse_datetime
        real(rk)       :: duration_s  = 0.0_rk   ! Duration in seconds
        real(rk)       :: t_start, t_end         ! Start and end time (in seconds since simulation start)
        ! Location
        real(rk) :: depth = 0._rk
        real(rk) :: min_depth, max_depth
        real(rk) :: thickness = 0._rk
        real(rk) :: thickness_eff = 0._rk      ! Effective thickness
        ! Indices of layers to apply the pulse
        integer  :: depth_idx = 0              
        integer  :: shallow_idx, deep_idx

        character(len=:), allocatable :: mode  ! "rate" or "amount"
        real(rk) :: rate   = 0._rk             ! mol m-3 s-1
        real(rk) :: amount = 0._rk             ! mol m-2 over duration 
        real(rk) :: amount_rate_vol = 0.0_rk   ! mol m-3 s-1
    end type

contains

    subroutine tracer_pulse_init(cfg, grid, FabmMod, sim_startdate, sim_enddate, cal, TP)
        type(ConfigParams),      intent(in)  :: cfg
        type(VerticalGrid),      intent(in)  :: grid
        class(type_fabm_model), pointer,  intent(in)  :: FabmMod
        type(DateTime),          intent(in)  :: sim_startdate, sim_enddate
        type(CFCalendar),        intent(in)  :: cal
        type(TracerPulse), intent(out) :: TP

        character(len=:), allocatable :: tracer, mode
        character(len=:), allocatable :: pulse_date_str, duration_str
        integer     :: idx, i, k_deep, k_shallow
        integer(lk) :: duration_sec, start_pulse_sec, sim_duration
        real(rk)    :: acc, dz_lo, dz_hi

        logical :: ok
        character(len=256) :: msg


        ! --- Location
        if (.not. cfg%is_set('location.depth')) call fatal('tracer_pulse:init', '"location.depth" missing.')
        TP%depth = cfg%get_param_num('location.depth', required=.true., min=0.0_rk)
        TP%thickness = cfg%get_param_num('location.thickness', default=0._rk, min=0.0_rk)

        if (TP%depth > grid%depth) call fatal('tracer_pulse:init', '"location.depth" exceeds location depth.')
        if (TP%thickness > grid%depth) call fatal('tracer_pulse:init', '"location.thickness" exceeds location depth.')

        TP%depth_idx = find_nearest_index(grid%z, TP%depth, prefer_min=.true.)
        TP%depth = grid%z(TP%depth_idx)

        ! Finding effective thickness to spread the pulse and corresponding layer indices
        if (TP%thickness > 0.0_rk) then        
            k_deep = TP%depth_idx
            k_shallow = TP%depth_idx
            acc  = grid%dz(TP%depth_idx)

            do while (acc < TP%thickness .and. (k_deep > 1 .or. k_shallow < grid%nz))
                ! Candidate expansion thicknesses 
                dz_lo = -1.0_rk
                dz_hi = -1.0_rk
                if (k_deep > 1)  dz_lo = grid%dz(k_deep-1)
                if (k_shallow < grid%nz) dz_hi = grid%dz(k_shallow+1)

                ! Expand towards the side that adds less overshoot
                if (dz_lo < 0.0_rk) then
                    k_shallow = k_shallow + 1
                    acc  = acc + dz_hi
                else if (dz_hi < 0.0_rk) then
                    k_deep = k_deep - 1
                    acc  = acc + dz_lo
                else
                    if (dz_lo <= dz_hi) then
                        k_deep = k_deep - 1
                        acc  = acc + dz_lo
                    else
                        k_shallow = k_shallow + 1
                        acc  = acc + dz_hi
                    end if
                end if
            end do

            TP%shallow_idx = k_shallow
            TP%deep_idx = k_deep

            TP%min_depth = grid%z_w(k_shallow)                 ! lower interface of deepest included layer
            TP%max_depth = grid%z_w(k_deep-1)                  ! upper interface of shallowest included layer
            TP%thickness_eff = TP%max_depth - TP%min_depth     ! Effective thickness

            if (abs(TP%thickness_eff - TP%thickness) > 1.0e-10_rk) then
                write(*,'(A)') 'WARNING [tracer_pulse]: requested thickness ('// &
                                trim(realtostr(TP%thickness,2))//'m) adjusted to match layer geometry.'
                write(*,'(A)') '         Effective thickness = '// &
                                trim(realtostr(TP%thickness_eff,2))//'m.'

            end if

        else 
            TP%thickness_eff = grid%dz(TP%depth_idx)
            TP%min_depth = grid%z_w(TP%depth_idx)
            TP%max_depth = grid%z_w(TP%depth_idx-1)
        end if

        !--- Tracer and mode
        tracer = cfg%get_param_str('events.tracer', required=.true., trim_value=.true.)
        mode   = to_lower(cfg%get_param_str('events.mode', default='rate', trim_value=.true.))

        if (mode /= 'rate' .and. mode /= 'amount') then
            call fatal('tracer_pulse:init', '"events.mode" must be "rate" or "amount".')
        end if
        TP%mode = mode
        TP%tracer_name = tracer

        ! ---- Find tracer index ----
        idx = find_tracer_index_by_name(FabmMod, tracer)
        if (idx < 1 .or. idx > size(FabmMod%interior_state_variables)) then
            write(*,'(A)') '[FATAL:tracer_pulse:init] Unknown tracer "'//trim(tracer)//'". Available interior tracers:'
            do i=1, size(FabmMod%interior_state_variables)
                write(*,'(A)') '  - '//trim(FabmMod%interior_state_variables(i)%name)
            end do
            stop 1
        end if
        TP%tracer_idx = idx  

        ! ---- require time fields ----
        if (.not. cfg%is_set('events.start'))    call fatal('tracer_pulse:init', '"events.start" missing.')
        if (.not. cfg%is_set('events.duration')) call fatal('tracer_pulse:init', '"events.duration" missing.')

        pulse_date_str = cfg%get_param_str('events.start', required=.true., trim_value=.true.)
        !Parse & validate dates according to the simulation calendar
        TP%pulse_datetime = datetime_from_string(trim(pulse_date_str), ok, msg, cal%kind)

        if (.not. is_datetime_before(TP%pulse_datetime, sim_enddate)) then
            if (.not. is_datetime_equal(TP%pulse_datetime, sim_enddate)) then
                write(*,'(A)') 'Invalid tracer pulse datetime: pulse occurs after the simulation end date.'
                stop 1
            end if
        end if
        if (.not. is_datetime_before(sim_startdate, TP%pulse_datetime) ) then
            if (.not. is_datetime_equal(sim_startdate, TP%pulse_datetime)) then 
                write(*,'(A)') 'Invalid tracer pulse datetime: pulse occurs before the simulation start date.'
                stop 1
            end if
        end if

        !-- Parse duration 
        duration_str = cfg%get_param_str('events.duration', default='1h', trim_value=.true.)
        call parse_interval_to_seconds(duration_str, duration_sec)   ! e.g., "1h" -> 3600
        if (duration_sec <=0_lk) then
           write(*,'(A)') 'events.duration is invalid; use positive forms like 900s, 15m, 1h, 1d'
           stop 1
        end if
        TP%duration_s = real(duration_sec, rk)

        start_pulse_sec = seconds_between_datetimes(cal, TP%pulse_datetime, sim_startdate)
        sim_duration    = seconds_between_datetimes(cal, sim_enddate, sim_startdate)
        TP%t_start = real(start_pulse_sec, rk)
        TP%t_end = TP%t_start + TP%duration_s
        if (TP%t_end > real(sim_duration,rk)) then
            write(*,'(A)') 'WARNING: events.duration extends beyond the end of simulation.'
        end if

        ! ---- magnitude requirements ----
        if (TP%mode == 'rate') then
            if (.not. cfg%is_set('events.rate')) call fatal('tracer_pulse:init', '"events.rate" required when mode=rate.')
            TP%rate = cfg%get_param_num('events.rate', required=.true., finite=.true., positive=.true.)
        else
            if (.not. cfg%is_set('events.amount')) call fatal('tracer_pulse:init', '"events.amount" required when mode=amount.')
            TP%amount = cfg%get_param_num('events.amount', required=.true., finite=.true., positive=.true.)
            ! Spread "amount [mol m-2]" over duration and thickness_eff
            TP%amount_rate_vol = (TP%amount / TP%duration_s) / TP%thickness_eff
        end if

    end subroutine tracer_pulse_init

    subroutine tracer_pulse_apply(model_time, dt, TP, tendency)        
        real(rk), intent(in)             :: model_time, dt
        type(TracerPulse), intent(in)    :: TP        
        real(rk),          intent(inout) :: tendency(:,:)   ! Tendency array for all tracers

        real(rk) :: t0, t1, overlap, frac
        real(rk) :: dC_over_dt
        integer  :: k, it

        t0 = model_time          ! Model time in seconds since the start of simulation.
        t1 = model_time + dt     ! Model time by the end of the current time-step
        overlap = max(0._rk, min(t1, TP%t_end) - max(t0, TP%t_start)) ! Overlap in seconds between the current time-step and the event interval
        if (overlap <= 0._rk) return

        ! Tracer index
        it = TP%tracer_idx

        frac = overlap / dt
        frac = max(0._rk, min(1._rk, frac))   ! Frac is between 0 and 1
        select case (TP%mode)
        case ('rate')
            ! Source active for overlap seconds: dC = rate * overlap
            ! frac is used to modulate whether the pulse is active for the whole substep
            ! If it’s active for half the substep, then apply half the rate over the whole substep.
            dC_over_dt = TP%rate * frac
        case ('amount')
            dC_over_dt = TP%amount_rate_vol * frac  ! mol m-3 s-1

        case default
            return  
        end select

        do k = TP%deep_idx, TP%shallow_idx
            tendency(k, it) = tendency(k, it) + dC_over_dt
        end do
    end subroutine tracer_pulse_apply


    subroutine tracer_pulse_clear(TP)
        type(TracerPulse), intent(inout) :: TP
        if (allocated(TP%tracer_name)) deallocate(TP%tracer_name)
        if (allocated(TP%mode))        deallocate(TP%mode)
        ! reset scalars/indices too
    end subroutine


    subroutine fatal(where, msg)
        character(len=*), intent(in) :: where, msg
        write(*,'(A,1X,A)') '[FATAL:'//trim(where)//']', trim(msg)
        stop 1
    end subroutine fatal

    pure integer function find_tracer_index_by_name(FabmMod, name) result(idx)
        class(type_fabm_model), pointer, intent(in) :: FabmMod
        character(len=*),       intent(in) :: name

        integer :: i
        character(len=:), allocatable :: q

        idx = -1
        q = to_lower(trim(adjustl(name)))

        do i = 1, size(FabmMod%interior_state_variables)
            if (to_lower(trim(FabmMod%interior_state_variables(i)%name)) == q) then
                idx = i
                return
            end if
        end do
    end function find_tracer_index_by_name

end module tracer_pulse
