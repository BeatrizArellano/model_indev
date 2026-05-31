module initial_profiles 
    use array_utils,       only: is_monotonic_increasing, is_monotonic_decreasing 
    use precision_utils,   only: is_equal_array
    use find_utils,        only: has_name
    use grids,             only: VerticalGrid, default_profile_depth_tol
    use path_utils,        only: path_has_extension
    use precision_types,   only: rk
    use state_loader,      only: StateData, build_state_from_tidy_table
    use str_utils,         only: inttostr, to_upper
    use text_parser,       only: TextTable, parse_text_table
    use variable_registry, only: VarMetadata

    implicit none
    private

    public :: set_initial_profiles

    integer, parameter  :: STATE_FORMAT_UNKNOWN = 0
    integer, parameter  :: STATE_FORMAT_TIDY    = 1
    integer, parameter  :: STATE_FORMAT_WIDE    = 2

    integer, parameter  :: NAME_LEN = 256
    integer, parameter  :: ERR_LEN  = 1024

contains 

    subroutine set_initial_profiles(filename, vars, grid, ok, errmsg)
        character(*),       intent(in)    :: filename
        type(VarMetadata),  intent(inout) :: vars(:)
        type(VerticalGrid), intent(in)    :: grid
        logical, optional,  intent(out)   :: ok
        character(len=*), optional, intent(out) :: errmsg

        type(StateData) :: state
        logical :: lok
        character(len=ERR_LEN) :: lmsg

        if (present(ok)) ok = .false.
        if (present(errmsg)) errmsg = ''

        call load_profile_data(filename, state, lok, lmsg)
        if (.not. lok) then
            call state%free()
            if (present(errmsg)) errmsg = trim(lmsg)
            if (.not. present(ok)) error stop trim(lmsg)
            return
        end if

        call apply_initial_profiles(state, vars, grid, lok, lmsg)
        if (.not. lok) then
            call state%free()
            if (present(errmsg)) errmsg = trim(lmsg)
            if (.not. present(ok)) error stop trim(lmsg)
            return
        end if

        call state%free()
        if (present(ok)) ok = .true.
    end subroutine set_initial_profiles



    !=========================================================
    !       Internal
    !=========================================================
    subroutine load_profile_data(filename, state, ok, errmsg)
        character(*),    intent(in)    :: filename
        type(StateData), intent(inout) :: state
        logical,         intent(out)   :: ok
        character(*),    intent(out)   :: errmsg

        type(TextTable) :: table
        integer  :: i_depth, i_name, i_value, fmt
        logical  :: exists, supported

        logical         :: lok
        character(len=ERR_LEN) :: lmsg

        exists = .false.
        supported = .false.
        ok = .false.
        errmsg = ''

        call state%free()

        if (len_trim(filename) == 0) then
            errmsg = 'Initial profile file path is empty.'
            return
        end if

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            errmsg = 'Initial profile file does not exist: '//trim(filename)
            return
        end if

        ! Verify extension
        supported = path_has_extension(filename, 'csv') .or. &
                    path_has_extension(filename, 'txt') .or. &
                    path_has_extension(filename, 'dat')

        if (.not. supported) then
            errmsg = 'Unsupported state/profile file extension: '//trim(filename)// &
                     '. Expected .csv, .txt, or .dat. NetCDF restart files are not supported yet.'
            return
        end if
        ! Parse text file into a table
        call parse_text_table(trim(filename), table)

        ! Figure out if it is wide or tidy format
        i_name  = table%find_column('name')
        i_depth = table%find_column('depth')
        i_value = table%find_column('value')

        fmt = STATE_FORMAT_UNKNOWN

        if (i_name > 0 .and. i_depth > 0 .and. i_value > 0) then
            fmt = STATE_FORMAT_TIDY
        else if (i_depth > 0 .and. table%ncol > 1) then
            fmt = STATE_FORMAT_WIDE
        end if
        if (fmt == STATE_FORMAT_UNKNOWN) then
            errmsg = 'Unsupported state/profile file format: '//trim(filename)// &
                     '. Profile file must be either tidy (name,depth,value) or wide (depth,var1,var2,...).'
            call table%free()
            return
        end if

        select case (fmt)
        case (STATE_FORMAT_TIDY)
            call build_state_from_tidy_table(table, state, lok, lmsg)
        case (STATE_FORMAT_WIDE)
            call build_profile_from_wide_table(table, state, lok, lmsg)
        end select

        if (.not. lok) then
            errmsg = trim(lmsg)
            call table%free()
            call state%free()
            return
        end if

        state%source_file = trim(filename)

        call table%free()
        ok = .true.
    end subroutine load_profile_data

    subroutine build_profile_from_wide_table(table, state, ok, errmsg)
        type(TextTable), intent(in)    :: table
        type(StateData), intent(inout) :: state
        logical,         intent(out)   :: ok
        character(*),    intent(out)   :: errmsg

        integer :: i_depth, jcol, ivar, irow, ios
        real(rk) :: z, val
        character(len=NAME_LEN) :: vname

        ok = .false.
        errmsg = ''

        i_depth = table%find_column('depth')

        if (i_depth <= 0) then
            errmsg = 'Wide state/profile file is missing required column: depth'
            return
        end if

        if (table%ncol <= 1) then
            errmsg = 'Wide state/profile file must contain depth and at least one variable column.'
            return
        end if

        if (table%nrow <= 0) then
            errmsg = 'Wide state/profile file contains no data rows.'
            return
        end if

        state%nvars = table%ncol - 1
        allocate(state%vars(state%nvars))

        ivar = 0
        do jcol = 1, table%ncol
            if (jcol == i_depth) cycle

            vname = trim(adjustl(table%header(jcol)))

            if (len_trim(vname) == 0) then
                errmsg = 'Empty variable name in wide state/profile file header.'
                return
            end if

            if (jcol > 1) then
                if (has_name(table%header(:jcol-1), vname)) then
                    errmsg = 'Duplicate variable name in wide state/profile file header: '//trim(vname)
                    return
                end if
            end if

            ivar = ivar + 1

            state%vars(ivar)%name = trim(vname)
            state%vars(ivar)%vert_coord = 'centre'      ! Setting to centre because only interior variables are expected here
            state%vars(ivar)%has_vert_coord = .true.
            state%vars(ivar)%nvals = table%nrow

            allocate(state%vars(ivar)%depth(table%nrow))
            allocate(state%vars(ivar)%values(table%nrow))

            do irow = 1, table%nrow
                read(table%values(irow, i_depth), *, iostat=ios) z
                if (ios /= 0) then
                    errmsg = 'Invalid depth value at line '// &
                            trim(inttostr(table%line_numbers(irow)))// &
                            ': '//trim(table%values(irow, i_depth))
                    return
                end if

                if (z < 0.0_rk) then
                    errmsg = 'Negative depth in wide state/profile file at line '// &
                            trim(inttostr(table%line_numbers(irow)))//'.'
                    return
                end if

                read(table%values(irow, jcol), *, iostat=ios) val
                if (ios /= 0) then
                    errmsg = 'Invalid state/profile value for variable '//trim(vname)// &
                            ' at line '//trim(inttostr(table%line_numbers(irow)))// &
                            ': '//trim(table%values(irow, jcol))
                    return
                end if

                state%vars(ivar)%depth(irow)  = z
                state%vars(ivar)%values(irow) = val
            end do
        end do

        ok = .true.
    end subroutine build_profile_from_wide_table

    subroutine apply_initial_profiles(state, vars, grid, ok, errmsg)
        type(StateData),    intent(in)    :: state
        type(VarMetadata),  intent(inout) :: vars(:)
        type(VerticalGrid), intent(in)    :: grid
        logical,            intent(out)   :: ok
        character(*),       intent(out)   :: errmsg

        real(rk), allocatable :: tol(:)
        integer  :: is, ivar, n_expected
        real(rk), allocatable :: depth_work(:), value_work(:)

        ok = .false.
        errmsg = ''

        n_expected = grid%nz

        if (state%nvars <= 0) then
            errmsg = 'Initial profile file contains no variables.'
            return
        end if

        if (grid%nz <= 0) then
            errmsg = 'Cannot set initial profiles on an empty vertical grid.'
            return
        end if

        call default_profile_depth_tol(grid, tol)

        do is = 1, state%nvars
            ivar = find_var_by_name(vars, state%vars(is)%name)

            if (ivar <= 0) then
                errmsg = 'Initial profile variable is not registered in interior variables: '// &
                        trim(state%vars(is)%name)
                return
            end if

            if (.not. vars(ivar)%state_var) then
                errmsg = 'Initial profile variable is not a state variable: '// &
                        trim(state%vars(is)%name)
                return
            end if

            if (trim(vars(ivar)%vert_coord) /= 'centre') then
                errmsg = 'Initial profile variable is not an interior variable: '// &
                          trim(state%vars(is)%name)// &
                          '. Domain is: '//trim(vars(ivar)%vert_coord)
                return
            end if

            if (state%vars(is)%has_vert_coord) then
                if (trim(state%vars(is)%vert_coord) /= 'centre') then
                    errmsg = 'Initial profile variable must use centre vertical coordinate: '// &
                              trim(state%vars(is)%name)// &
                             '. File vertical coordinate is: '//trim(state%vars(is)%vert_coord)
                    return
                end if
            end if

            if (.not. associated(vars(ivar)%data_1d)) then
                errmsg = 'Initial profile variable has no associated 1D storage: '// &
                        trim(state%vars(is)%name)
                return
            end if

            if (state%vars(is)%nvals /= n_expected) then
                errmsg = 'Initial profile variable '//trim(state%vars(is)%name)// &
                        ' has wrong number of values. Expected '//trim(inttostr(n_expected))// &
                        ', got '//trim(inttostr(state%vars(is)%nvals))//'.'
                return
            end if

            if (allocated(depth_work)) deallocate(depth_work)
            if (allocated(value_work)) deallocate(value_work)
            allocate(depth_work(n_expected), value_work(n_expected))

            if (is_monotonic_decreasing(state%vars(is)%depth, 1.0e-6_rk)) then
                depth_work = state%vars(is)%depth
                value_work = state%vars(is)%values

            else if (is_monotonic_increasing(state%vars(is)%depth, 1.0e-6_rk)) then
                depth_work = state%vars(is)%depth(n_expected:1:-1)
                value_work = state%vars(is)%values(n_expected:1:-1)

            else
                errmsg = 'Depth coordinate for initial profile variable '// &
                         trim(state%vars(is)%name)//' is not monotonic.'
                return
            end if

            if (.not. all(is_equal_array(depth_work, grid%z, abs_=tol))) then
                errmsg = 'Depth coordinate for initial profile variable '// &
                          trim(state%vars(is)%name)//' does not match model layer-centre grid.'
                if (allocated(tol)) deallocate(tol)
                return
            end if

            vars(ivar)%data_1d(:) = value_work(:)

        end do

        if (allocated(depth_work)) deallocate(depth_work)
        if (allocated(value_work)) deallocate(value_work)
        if (allocated(tol)) deallocate(tol)

        ok = .true.
    end subroutine apply_initial_profiles

    integer function find_var_by_name(vars, name) result(idx)
        type(VarMetadata), intent(in) :: vars(:)
        character(*),      intent(in) :: name

        integer :: i

        idx = -1
        do i = 1, size(vars)
            if (to_upper(trim(adjustl(vars(i)%name))) == to_upper(trim(adjustl(name)))) then
                idx = i
                return
            end if
        end do
    end function find_var_by_name    

end module initial_profiles  


    