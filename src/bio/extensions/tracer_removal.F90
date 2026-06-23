! Notes:
!   - This event removes an existing water-column tracer by directly modifying state.
!   - Application regions are represented by complete model layers.
!   - thickness_eff is the total thickness of the complete model layers selected
!     for application.
!   - Areal removals are converted internally to concentration changes by
!     dividing by thickness_eff.
module tracer_removal
    use bio_types,        only: BioEnv
    use event_types,      only: Event, EventData
    use event_utils,      only: find_tracer_index_by_name, resolve_application_layers, fatal
    use precision_types,  only: rk
    use read_config_yaml, only: ConfigParams
    use str_utils,        only: to_lower

    implicit none
    private

    public :: TracerRemoval
    public :: tracer_removal_prepare, tracer_removal_apply, tracer_removal_clear

    type, extends(EventData) :: TracerRemoval

        ! Tracer information
        integer :: tracer_idx = -1
        character(len=:), allocatable :: tracer_name

        ! Location: surface, bottom, or depth
        character(len=:), allocatable :: location_kind
        real(rk) :: depth         = 0._rk
        real(rk) :: min_depth     = 0._rk
        real(rk) :: max_depth     = 0._rk
        real(rk) :: thickness     = 0._rk
        real(rk) :: thickness_eff = 0._rk
        integer  :: depth_idx     = 0
        integer  :: shallow_idx   = 0
        integer  :: deep_idx      = 0

        ! Removal mode:
        character(len=:), allocatable :: mode  ! "rate" or "amount"

        real(rk) :: rate_conc = 0._rk   ! concentration/time sink
        real(rk) :: rate_area = 0._rk   ! areal/time removal
        real(rk) :: amount    = 0._rk   ! total areal amount to remove over duration

        ! Diagnostics/bookkeeping in areal units
        real(rk) :: amount_initial   = 0._rk
        real(rk) :: amount_remaining = 0._rk
        real(rk) :: amount_removed   = 0._rk
    end type TracerRemoval

contains

    subroutine tracer_removal_prepare(evt, cfg, BE)
        type(Event),        intent(inout) :: evt
        type(ConfigParams), intent(in)    :: cfg
        type(BioEnv),       intent(in)    :: BE

        character(len=:), allocatable :: basekey, tracer, mode, location
        integer :: idx, i

        if (allocated(evt%details)) deallocate(evt%details)
        allocate(TracerRemoval :: evt%details)

        select type (TR => evt%details)
        type is (TracerRemoval)

            basekey = 'events.' // trim(evt%name)

            !--------------------------------------------------
            ! This event directly modifies the tracer state.
            ! It does not contribute to BE%tendency_int_evt.
            !--------------------------------------------------
            evt%changes_tendency       = .false.
            evt%changes_state_directly = .true.
            evt%is_instantaneous       = .false.

            if (evt%duration <= 0._rk) then
                call fatal('tracer_removal:prepare', &
                    'Event "'//trim(evt%name)//'" duration should be positive.')
            end if

            !--------------------------------------------------
            ! Location: water column only
            !--------------------------------------------------
            location = to_lower(cfg%get_param_str(basekey//'.location', default='depth', trim_value=.true.))

            TR%location_kind = location
            TR%thickness = cfg%get_param_num(basekey//'.thickness', default=0._rk, min=0.0_rk)

            call resolve_application_layers(basekey, location, cfg, BE%wat_grid,                    &
                                            TR%depth, TR%depth_idx, TR%thickness, TR%thickness_eff, &
                                            TR%shallow_idx, TR%deep_idx, TR%min_depth, TR%max_depth)

            !--------------------------------------------------
            ! Tracer and mode
            !--------------------------------------------------
            tracer = cfg%get_param_str(basekey//'.tracer', required=.true., trim_value=.true.)
            mode = to_lower(cfg%get_param_str(basekey//'.mode', default='rate', trim_value=.true.))

            if (mode /= 'rate' .and. mode /= 'amount') then
                call fatal('tracer_removal:prepare', '"'//basekey//'.mode" must be "rate" or "amount".')
            end if

            TR%mode        = mode
            TR%tracer_name = tracer

            idx = find_tracer_index_by_name(BE%model, tracer)

            if (idx < 1 .or. idx > size(BE%model%interior_state_variables)) then
                write(*,'(A)') '[FATAL:tracer_removal:prepare] Unknown tracer "'//trim(tracer)//'". Available interior tracers:'
                do i = 1, size(BE%model%interior_state_variables)
                    write(*,'(A)') '  - '//trim(BE%model%interior_state_variables(i)%name)
                end do
                stop 1
            end if

            TR%tracer_idx = idx

            !--------------------------------------------------
            ! Magnitude
            !--------------------------------------------------
            select case (TR%mode)
            case ('rate')
                if (.not. cfg%is_set(basekey//'.rate')) then
                    call fatal('tracer_removal:prepare', '"'//trim(basekey)//'.rate" required when mode="rate".')
                end if

                TR%rate_area = cfg%get_param_num(basekey//'.rate', required=.true., finite=.true., positive=.true.)
                TR%rate_conc = TR%rate_area / TR%thickness_eff

                TR%amount_initial   = 0._rk
                TR%amount_remaining = 0._rk
                TR%amount_removed   = 0._rk

            case ('amount')
                if (.not. cfg%is_set(basekey//'.amount')) then
                    call fatal('tracer_removal:prepare', '"'//trim(basekey)//'.amount" required when mode="amount".')
                end if

                TR%amount = cfg%get_param_num(basekey//'.amount', required=.true., finite=.true., positive=.true.)

                TR%amount_initial   = TR%amount
                TR%amount_remaining = TR%amount
                TR%amount_removed   = 0._rk
                ! For mode="amount", the instantaneous removal rate is not fixed here.
            end select            

        end select

    end subroutine tracer_removal_prepare


    subroutine tracer_removal_apply(evt, BE, model_time, dt)
        type(Event),  intent(inout) :: evt
        type(BioEnv), intent(inout) :: BE
        real(rk),     intent(in)    :: model_time
        real(rk),     intent(in)    :: dt

        integer  :: k, it
        real(rk) :: requested_dC, available_dC, removed_dC
        real(rk) :: requested_area_step, removed_area_step
        real(rk) :: time_remaining, dt_eff, rate_area_now
        real(rk) :: lower_bound

        if (.not. evt%changes_state_directly) return
        if (model_time + dt <= evt%t_start) return
        if (model_time >= evt%t_end) return

        select type (TR => evt%details)
        type is (TracerRemoval)

            it = TR%tracer_idx

            time_remaining = evt%t_end - model_time
            if (time_remaining <= 0._rk) return

            dt_eff = min(dt, time_remaining)

            select case (TR%mode)

            case ('rate')
                ! Fixed removal rate per squared meter:
                !   rate_area [tracer units m-2 s-1]
                rate_area_now = TR%rate_area

            case ('amount')
                ! Adaptive areal removal rate:
                ! remove the remaining requested amount over the remaining event time.
                if (TR%amount_remaining <= 0._rk) return

                rate_area_now = TR%amount_remaining / time_remaining

            case default
                call fatal('tracer_removal:apply', &
                    'Unknown removal mode "'//trim(TR%mode)//'".')

            end select

            requested_area_step = rate_area_now * dt_eff

            if (TR%mode == 'amount') then
                requested_area_step = min(requested_area_step, TR%amount_remaining)
            end if

            requested_dC = requested_area_step / TR%thickness_eff
            if (requested_dC <= 0._rk) return

            lower_bound = BE%model%interior_state_variables(it)%minimum
            removed_area_step = 0._rk

            do k = TR%deep_idx, TR%shallow_idx
                available_dC = max(0._rk, BE%BS%interior_state(k, it) - lower_bound)
                removed_dC   = min(requested_dC, available_dC)

                BE%BS%interior_state(k, it) = BE%BS%interior_state(k, it) - removed_dC

                ! Convert concentration removed from this layer to areal inventory.
                removed_area_step = removed_area_step + removed_dC * BE%wat_grid%dz(k)
            end do

            TR%amount_removed = TR%amount_removed + removed_area_step

            if (TR%mode == 'amount') then
                TR%amount_remaining = max(0._rk, TR%amount_remaining - removed_area_step)
            end if

        class default
            call fatal('tracer_removal:apply', 'Invalid event details type.')

        end select
    end subroutine tracer_removal_apply


    subroutine tracer_removal_clear(TR)
        type(TracerRemoval), intent(inout) :: TR

        if (allocated(TR%tracer_name))   deallocate(TR%tracer_name)
        if (allocated(TR%mode))          deallocate(TR%mode)
        if (allocated(TR%location_kind)) deallocate(TR%location_kind)
    end subroutine tracer_removal_clear

end module tracer_removal