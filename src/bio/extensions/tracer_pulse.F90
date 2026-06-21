module tracer_pulse
   use bio_types,        only: BioEnv
   use event_types,      only: Event, EventData
   use fabm,             only: type_fabm_model
   use find_utils,       only: find_nearest_index
   use precision_types,  only: rk
   use read_config_yaml, only: ConfigParams
   use str_utils,        only: to_lower

   implicit none
   private
   public :: tracer_pulse_prepare, tracer_pulse_apply, tracer_pulse_clear
   public :: TracerPulse

    type, extends(EventData) :: TracerPulse
        ! Tracer
        integer :: tracer_idx = -1
        character(len=:), allocatable :: tracer_name
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

    subroutine tracer_pulse_prepare(evt, cfg, BE)
        type(Event),        intent(inout) :: evt
        type(ConfigParams), intent(in)    :: cfg
        type(BioEnv),       intent(in)    :: BE

        character(len=:), allocatable :: basekey, tracer, mode
        integer  :: idx, i, k_deep, k_shallow
        real(rk) :: acc, dz_lo, dz_hi

        if (allocated(evt%details)) deallocate(evt%details)
        allocate(TracerPulse :: evt%details)

        select type (TP => evt%details)
        type is (TracerPulse)

            basekey = 'events.' // trim(evt%name)

            !--------------------------------------------------
            ! This event changes water-column tracer tendencies
            !--------------------------------------------------
            evt%changes_tendency = .true.
            ! It is applied over a time interval
            evt%is_instantaneous = .false.   
            
            if (evt%duration <= 0._rk) then
                call fatal('tracer_pulse:prepare', &
                        'Event "'//trim(evt%name)//'" has non-positive duration.')
            end if

            !--------------------------------------------------
            ! Location: water column only
            !--------------------------------------------------
            if (.not. cfg%is_set(basekey//'.depth')) then
                call fatal('tracer_pulse:prepare', '"'//basekey//'.depth" missing.')
            end if

            TP%depth = cfg%get_param_num(basekey//'.depth', required=.true., min=0.0_rk)
            TP%thickness = cfg%get_param_num(basekey//'.thickness', default=0._rk, min=0.0_rk)

            if (TP%depth > BE%wat_grid%depth) then
                call fatal('tracer_pulse:prepare', '"'//basekey//'.depth" exceeds water-column depth.')
            end if

            if (TP%thickness > BE%wat_grid%depth) then
                call fatal('tracer_pulse:prepare', '"'//basekey//'.thickness" exceeds water-column depth.')
            end if

            TP%depth_idx = find_nearest_index(BE%wat_grid%z, TP%depth, prefer_min=.true.)
            TP%depth = BE%wat_grid%z(TP%depth_idx)

            ! Finding effective thickness to spread the pulse and corresponding layer indices
            if (TP%thickness > 0.0_rk) then
                k_deep = TP%depth_idx
                k_shallow = TP%depth_idx
                acc = BE%wat_grid%dz(TP%depth_idx)

                do while (acc < TP%thickness .and. (k_deep > 1 .or. k_shallow < BE%wat_grid%nz))
                    ! Candidate expansion thicknesses 
                    dz_lo = -1.0_rk
                    dz_hi = -1.0_rk

                    if (k_deep > 1) dz_lo = BE%wat_grid%dz(k_deep-1)
                    if (k_shallow < BE%wat_grid%nz) dz_hi = BE%wat_grid%dz(k_shallow+1)

                    if (dz_lo < 0.0_rk) then
                        k_shallow = k_shallow + 1
                        acc = acc + dz_hi
                    else if (dz_hi < 0.0_rk) then
                        k_deep = k_deep - 1
                        acc = acc + dz_lo
                    else if (dz_lo <= dz_hi) then
                        k_deep = k_deep - 1
                        acc = acc + dz_lo
                    else
                        k_shallow = k_shallow + 1
                        acc = acc + dz_hi
                    end if
                end do

                TP%shallow_idx = k_shallow
                TP%deep_idx = k_deep

                TP%min_depth = BE%wat_grid%z_w(k_shallow)
                TP%max_depth = BE%wat_grid%z_w(k_deep-1)
                TP%thickness_eff = TP%max_depth - TP%min_depth

            else
                TP%shallow_idx   = TP%depth_idx
                TP%deep_idx      = TP%depth_idx
                TP%thickness_eff = BE%wat_grid%dz(TP%depth_idx)
                TP%min_depth     = BE%wat_grid%z_w(TP%depth_idx)
                TP%max_depth     = BE%wat_grid%z_w(TP%depth_idx-1)
            end if

            if (TP%thickness_eff <= 0._rk) then
                call fatal('tracer_pulse:prepare', &
                    'Computed non-positive effective thickness for event "'//trim(evt%name)//'".')
            end if

            !--------------------------------------------------
            ! Tracer and mode
            !--------------------------------------------------
            tracer = cfg%get_param_str(basekey//'.tracer', required=.true., trim_value=.true.)
            mode = to_lower(cfg%get_param_str(basekey//'.mode', default='rate', trim_value=.true.))

            if (mode /= 'rate' .and. mode /= 'amount') then
                call fatal('tracer_pulse:prepare', '"'//basekey//'.mode" must be "rate" or "amount".')
            end if

            TP%mode = mode
            TP%tracer_name = tracer

            ! ---- Find tracer index ----
            idx = find_tracer_index_by_name(BE%model, tracer)

            if (idx < 1 .or. idx > size(BE%model%interior_state_variables)) then
                write(*,'(A)') '[FATAL:tracer_pulse:prepare] Unknown tracer "'//trim(tracer)//'". Available interior tracers:'
                do i = 1, size(BE%model%interior_state_variables)
                    write(*,'(A)') '  - '//trim(BE%model%interior_state_variables(i)%name)
                end do
                stop 1
            end if

            TP%tracer_idx = idx

            !--------------------------------------------------
            ! Magnitude
            !--------------------------------------------------
            if (TP%mode == 'rate') then
                if (.not. cfg%is_set(basekey//'.rate')) call fatal('tracer_pulse:prepare', '"rate" required when mode=rate.')
                TP%rate = cfg%get_param_num(basekey//'.rate', required=.true., finite=.true., positive=.true.)

            else
                if (.not. cfg%is_set(basekey//'.amount')) call fatal('tracer_pulse:prepare', '"amount" required when mode=amount.')
                TP%amount = cfg%get_param_num(basekey//'.amount', required=.true., finite=.true., positive=.true.)
                TP%amount_rate_vol = (TP%amount / evt%duration) / TP%thickness_eff
            end if
        end select

    end subroutine tracer_pulse_prepare

    subroutine tracer_pulse_apply(evt, BE)

        type(Event),  intent(in)    :: evt
        type(BioEnv), intent(inout) :: BE

        real(rk) :: dC_over_dt
        integer  :: k, it

        if (.not. evt%changes_tendency) return
        ! This routine assumes that the accepted substep does not cross
        ! evt%t_start or evt%t_end. Event-boundary timestep limiting is
        ! handled by EventManager / bio_main.
        ! Event tendencies are only active on the half-open interval:
        ! [t_start, t_end). 

        select type (TP => evt%details)
        type is (TracerPulse)

            it = TP%tracer_idx
            if (it < 1) call fatal('tracer_pulse:apply', &
                'Tracer index was not initialized for event "'//trim(evt%name)//'".')

            select case (TP%mode)
            case ('rate')
                ! Volumetric source rate [tracer units m-3 s-1].
                dC_over_dt = TP%rate 

            case ('amount')
                ! During prepare, the areal amount [tracer units m-2] was
                ! converted to a constant volumetric rate over the event
                ! duration and effective application thickness.
                dC_over_dt = TP%amount_rate_vol 
            end select

            ! Add to the existing tendencies from other events
            do k = TP%deep_idx, TP%shallow_idx
                BE%tendency_int_evt(k, it) = BE%tendency_int_evt(k, it) + dC_over_dt
            end do
        end select

    end subroutine tracer_pulse_apply  


    subroutine tracer_pulse_clear(TP)
        type(TracerPulse), intent(inout) :: TP
        if (allocated(TP%tracer_name)) deallocate(TP%tracer_name)
        if (allocated(TP%mode))        deallocate(TP%mode)
        ! reset scalars/indices too
    end subroutine

    !======================================================================
    !    Internal
    !======================================================================


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
