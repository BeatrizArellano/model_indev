! Notes:
!   - Application regions are represented by complete model layers.
!   - thickness_eff is the total thickness of the complete model layers
!     selected for application. It is always <= the requested thickness.
module tracer_pulse
   use bio_types,        only: BioEnv
   use event_types,      only: Event, EventData
   use event_utils,      only: find_tracer_index_by_name, resolve_application_layers, fatal
   use fabm,             only: type_fabm_model
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
        character(len=:), allocatable :: location_kind  ! surface, bottom, depth
        ! Location
        real(rk) :: depth = 0._rk
        real(rk) :: min_depth, max_depth
        real(rk) :: thickness = 0._rk
        real(rk) :: thickness_eff = 0._rk      ! Effective thickness
        ! Indices of layers to apply the pulse
        integer  :: depth_idx = 0              
        integer  :: shallow_idx, deep_idx

        character(len=:), allocatable :: mode  ! "rate" or "amount"
        real(rk) :: rate_area = 0._rk             ! tracer units m s-1 (e.g. mmol m-2 s-1)
        real(rk) :: amount    = 0._rk             ! tracer units m-2
        real(rk) :: rate_vol  = 0.0_rk            ! tracer concentration units s-1, e.g. mmol m-3 s-1
    end type

contains

    subroutine tracer_pulse_prepare(evt, cfg, BE)
        type(Event),        intent(inout) :: evt
        type(ConfigParams), intent(in)    :: cfg
        type(BioEnv),       intent(in)    :: BE

        character(len=:), allocatable :: basekey, tracer, mode, location
        integer  :: idx, i

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
            evt%is_instantaneous       = .false.  
            ! This event is applied via a tendency and does not modify state directly
            evt%changes_state_directly = .false.  
            
            if (evt%duration <= 0._rk) then
                call fatal('tracer_pulse:prepare', &
                        'Event "'//trim(evt%name)//'" has non-positive duration.')
            end if

            !--------------------------------------------------
            ! Location: water column only
            !--------------------------------------------------
            location = to_lower(cfg%get_param_str(basekey//'.location', default='depth', trim_value=.true.))

            TP%location_kind = location
            TP%thickness = cfg%get_param_num(basekey//'.thickness', default=0._rk, min=0.0_rk)

            call resolve_application_layers(basekey, location, cfg, BE%wat_grid,                    &
                                            TP%depth, TP%depth_idx, TP%thickness, TP%thickness_eff, &
                                            TP%shallow_idx, TP%deep_idx, TP%min_depth, TP%max_depth)

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
                if (.not. cfg%is_set(basekey//'.rate')) then
                    call fatal('tracer_pulse:prepare', '"rate" required when mode=rate.')
                end if

                TP%rate_area = cfg%get_param_num(basekey//'.rate', required=.true., finite=.true., positive=.true.)
                TP%rate_vol = TP%rate_area / TP%thickness_eff

            else
                if (.not. cfg%is_set(basekey//'.amount')) then
                    call fatal('tracer_pulse:prepare', '"amount" required when mode=amount.')
                end if

                TP%amount = cfg%get_param_num(basekey//'.amount', required=.true., finite=.true., positive=.true.)
                TP%rate_vol = (TP%amount / evt%duration) / TP%thickness_eff
            end if

        end select

    end subroutine tracer_pulse_prepare

    subroutine tracer_pulse_apply(evt, BE)

        type(Event),  intent(in)    :: evt
        type(BioEnv), intent(inout) :: BE

        integer  :: k, it

        if (.not. evt%changes_tendency) return
        ! Event tendencies are only active on the half-open interval:
        ! [t_start, t_end). 

        select type (TP => evt%details)
        type is (TracerPulse)

            it = TP%tracer_idx

            ! Add to the existing tendencies from other events
            do k = TP%deep_idx, TP%shallow_idx
                BE%tendency_int_evt(k, it) = BE%tendency_int_evt(k, it) + TP%rate_vol
            end do
        end select

    end subroutine tracer_pulse_apply  


    subroutine tracer_pulse_clear(TP)
        type(TracerPulse), intent(inout) :: TP
        if (allocated(TP%tracer_name))   deallocate(TP%tracer_name)
        if (allocated(TP%mode))          deallocate(TP%mode)
        if (allocated(TP%location_kind)) deallocate(TP%location_kind)
    end subroutine

end module tracer_pulse
