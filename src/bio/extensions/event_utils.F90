module event_utils
    use fabm,             only: type_fabm_model
    use find_utils,       only: find_nearest_index    
    use grids,            only: VerticalGrid
    use precision_types,  only: rk
    use read_config_yaml, only: ConfigParams
    use str_utils,        only: to_lower

    implicit none
    private

    public :: find_tracer_index_by_name
    public :: resolve_application_layers
    public :: fatal

contains    

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
    
    !------------------------------------------------------------------------------------------------
    !     Finds the layer indices to apply the requested event
    ! The resolved application region always contains complete model layers.
    ! If requested thickness is thinner than one layer, the event is applied
    ! over that whole layer and diluted by thickness_eff.
    ! The effective thickness is the closest small number to thickness determined by the grid layers.
    !-------------------------------------------------------------------------------------------------
    subroutine resolve_application_layers(basekey, location, cfg, grid, depth, depth_idx, thickness, &
                                          thickness_eff, shallow_idx, deep_idx, min_depth, max_depth)
        character(len=*),    intent(in)  :: basekey, location
        type(ConfigParams),  intent(in)  :: cfg
        type(VerticalGrid),  intent(in)  :: grid
        real(rk),            intent(in)  :: thickness
        real(rk),            intent(out) :: depth, thickness_eff, min_depth, max_depth
        integer,             intent(out) :: depth_idx, shallow_idx, deep_idx

        integer  :: k_deep, k_shallow
        real(rk) :: acc
        logical :: can_deep, can_shallow
        character(len=:), allocatable :: loc

        loc = to_lower(trim(adjustl(location)))

        if (loc /= 'surface' .and. loc /= 'bottom' .and. loc /= 'depth') then
            call fatal('event_utils:resolve_application_layers', &
                '"'//trim(basekey)//'.location" must be "surface", "bottom", or "depth".')
        end if

        select case (loc)
            case ('surface')
                depth_idx = grid%nz         ! Surface index
                depth = grid%z(depth_idx)

            case ('bottom')                 ! Bottom index is always 1 for layer centres
                depth_idx = 1
                depth = grid%z(depth_idx)

            case ('depth')
                if (.not. cfg%is_set(basekey//'.depth')) then 
                    call fatal('event_utils:resolve_application_layers', &
                        '"'//trim(basekey)//'.depth" is required when location="depth".')
                end if

                depth = cfg%get_param_num(basekey//'.depth', required=.true., finite=.true., min=0._rk)
                if (depth >  grid%depth) then 
                    call fatal('event_utils:resolve_application_layers', &
                        '"'//trim(basekey)//'.depth" is deeper than water depth in this location.')
                end if 

                depth_idx = find_nearest_index(grid%z, depth, prefer_min=.true.)
                depth     = grid%z(depth_idx)    ! Use new depth
        end select

        if (depth_idx < 1 .or. depth_idx > grid%nz) then
            call fatal('event_utils:resolve_application_layers', &
                'Resolved depth index is outside the water-column grid for "'//trim(basekey)//'".')
        end if

        ! If no thickness is provided, use only one layer.
        if (thickness <= 0._rk) then
            deep_idx      = depth_idx
            shallow_idx   = depth_idx
            min_depth     = grid%z_w(depth_idx)
            max_depth     = grid%z_w(depth_idx-1)
            thickness_eff = grid%dz(depth_idx)
            return
        end if

        ! For surface/bottom locations, extend inward without exceeding the requested thickness.
        if (loc == 'surface') then
            shallow_idx = grid%nz
            deep_idx    = grid%nz
            acc = grid%dz(grid%nz)

            do while (deep_idx > 1)
                if (acc + grid%dz(deep_idx-1) > thickness) exit
                deep_idx = deep_idx - 1
                acc = acc + grid%dz(deep_idx)
            end do

        else if (loc == 'bottom') then
            deep_idx = 1
            shallow_idx = 1
            acc = grid%dz(1)

            do while (shallow_idx < grid%nz)
                if (acc + grid%dz(shallow_idx+1) > thickness) exit
                shallow_idx = shallow_idx + 1
                acc = acc + grid%dz(shallow_idx)
            end do
        else
            ! For location=depth, grow around the nearest layer without exceeding the requested thickness.
            deep_idx = depth_idx
            shallow_idx = depth_idx
            acc = grid%dz(depth_idx)

            do
                can_deep = .false.
                can_shallow = .false.

                if (deep_idx > 1) then
                    can_deep = acc + grid%dz(deep_idx-1) <= thickness
                end if

                if (shallow_idx < grid%nz) then
                    can_shallow = acc + grid%dz(shallow_idx+1) <= thickness
                end if

                if (.not. can_deep .and. .not. can_shallow) exit

                if (can_deep .and. can_shallow) then
                    if (grid%dz(deep_idx-1) <= grid%dz(shallow_idx+1)) then
                        deep_idx = deep_idx - 1
                        acc = acc + grid%dz(deep_idx)
                    else
                        shallow_idx = shallow_idx + 1
                        acc = acc + grid%dz(shallow_idx)
                    end if
                else if (can_deep) then
                    deep_idx = deep_idx - 1
                    acc = acc + grid%dz(deep_idx)
                else
                    shallow_idx = shallow_idx + 1
                    acc = acc + grid%dz(shallow_idx)
                end if
            end do
        end if

        k_deep = min(deep_idx, shallow_idx)
        k_shallow = max(deep_idx, shallow_idx)

        deep_idx = k_deep
        shallow_idx = k_shallow

        min_depth = grid%z_w(shallow_idx)
        max_depth = grid%z_w(deep_idx-1)
        thickness_eff = max_depth - min_depth

        if (thickness_eff <= 0._rk) then
            call fatal('event_utils:resolve_application_layers', &
                       'Resolved non-positive application thickness for "'//trim(basekey)//'".')
        end if

    end subroutine resolve_application_layers

    subroutine fatal(where, msg)
        character(len=*), intent(in) :: where, msg
        write(*,'(A,1X,A)') '[FATAL:'//trim(where)//']', trim(msg)
        stop 1
    end subroutine fatal

end module event_utils