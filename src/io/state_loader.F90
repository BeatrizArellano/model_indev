!=====================================================================
! state_loader
!---------------------------------------------------------------------
! Loader for initial/restart model state from tidy text files.
!
! Supported formats for now:
!   - .csv
!   - .txt / .dat  whitespace-delimited
!
! Required tidy columns:
!   name, depth, value
!=====================================================================
module state_loader
    use array_utils,       only: is_monotonic_increasing, is_monotonic_decreasing
    use grids,             only: VerticalGrid, default_profile_depth_tol
    use precision_types,   only: rk
    use precision_utils,   only: is_equal_array
    use path_utils,        only: path_has_extension
    use str_utils,         only: inttostr, to_upper
    use text_parser,       only: TextTable, parse_text_table
    use variable_registry, only: VarMetadata

    implicit none
    private

    public :: StateVariable
    public :: StateData
    public :: load_state_file, set_initial_state, build_state_from_tidy_table

    integer, parameter :: NAME_LEN = 256
    integer, parameter :: ERR_LEN  = 1024

    type :: StateVariable
        character(len=:), allocatable :: name
        character(len=16)             :: vert_coord = ''   ! centre/surface/bottom/none
        integer                       :: nvals = 0
        real(rk), allocatable         :: depth(:)
        real(rk), allocatable         :: values(:)
        logical :: has_vert_coord = .false.
    end type StateVariable

    type :: StateData
        character(len=:), allocatable :: source_file
        type(StateVariable), allocatable :: vars(:)
        integer :: nvars = 0        
    contains
        procedure :: free => state_data_free
        procedure :: find => state_data_find
    end type StateData

contains

    subroutine load_state_file(filename, state, ok, errmsg)
        character(*),    intent(in)    :: filename
        type(StateData), intent(inout) :: state
        logical,         intent(out)   :: ok
        character(*),    intent(out)   :: errmsg

        type(TextTable) :: table
        logical         :: lok
        character(len=ERR_LEN) :: lmsg

        ok = .false.
        errmsg = ''

        call state%free()

        call validate_state_path(filename, lok, lmsg)
        if (.not. lok) then
            errmsg = trim(lmsg)
            return
        end if        

        ! Parse the file and load it into a table.
        call parse_text_table(trim(filename), table)

        call validate_state_columns(table, lok, lmsg)
        if (.not. lok) then
            errmsg = trim(lmsg)
            call table%free()
            return
        end if

        call build_state_from_tidy_table(table, state, lok, lmsg)
        if (.not. lok) then
            errmsg = trim(lmsg)
            call table%free()
            call state%free()
            return
        end if
        state%source_file = trim(filename)

        call table%free()

        ok = .true.
    end subroutine load_state_file

    integer function state_data_find(self, name) result(idx)
        class(StateData), intent(in) :: self
        character(*),     intent(in) :: name

        integer :: i

        idx = 0

        if (.not. allocated(self%vars)) return

        do i = 1, self%nvars
            if (to_upper(trim(adjustl(self%vars(i)%name))) == &
                to_upper(trim(adjustl(name)))) then
                idx = i
                return
            end if
        end do
    end function state_data_find


    subroutine state_data_free(self)
        class(StateData), intent(inout) :: self

        if (allocated(self%vars)) deallocate(self%vars)
        if (allocated(self%source_file)) deallocate(self%source_file)

        self%nvars = 0
    end subroutine state_data_free


    subroutine set_initial_state(state, vars, grid, ok, errmsg)
        type(StateData),              intent(in)    :: state
        type(VarMetadata),            intent(inout) :: vars(:)
        type(VerticalGrid), optional, intent(in)    :: grid
        logical,            optional, intent(out)   :: ok
        character(len=*),   optional, intent(out)   :: errmsg

        real(rk), allocatable :: tol(:)
        real(rk) :: mono_tol
        integer  :: ivar, is, n_expected
        real(rk), allocatable :: depth_work(:), value_work(:)
        logical  :: depth_match, is_centre

        if (present(ok)) ok = .false.
        if (present(errmsg)) errmsg = ''

        mono_tol = 1.0e-6_rk
       
        do ivar = 1, size(vars)
            if (.not. vars(ivar)%state_var) cycle
            ! Find state variable by name
            is = state%find(vars(ivar)%name)

            if (is <= 0) then
                if (present(errmsg)) errmsg = 'Missing required state variable in initial state file: '//trim(vars(ivar)%name)
                if (present(ok)) then
                    ok = .false.
                    return
                else
                    error stop 'Missing required state variable in initial state file: '//trim(vars(ivar)%name)
                end if
            end if

            ! If the state file provided a vertical coordinate label, it must match the registry.
            if (len_trim(state%vars(is)%vert_coord) > 0) then
                if (to_upper(trim(state%vars(is)%vert_coord)) /= to_upper(trim(vars(ivar)%vert_coord))) then
                    if (present(errmsg)) errmsg = 'Vertical coordinate mismatch for variable '//trim(vars(ivar)%name)// &
                                                ': registry='//trim(vars(ivar)%vert_coord)// &
                                                ', file='//trim(state%vars(is)%vert_coord)
                    if (present(ok)) then
                        ok = .false.
                        return
                    else
                        error stop 'set_initial_state: vertical coordinate mismatch for variable '//trim(vars(ivar)%name)
                    end if
                end if
            end if

            select case (trim(vars(ivar)%vert_coord))
                case ('surface', 'bottom')
                    if (.not. associated(vars(ivar)%data_0d)) then
                        if (present(errmsg)) errmsg = 'State variable '//trim(vars(ivar)%name)//' has no associated scalar storage.'
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: unassociated scalar storage for '//trim(vars(ivar)%name)
                        end if
                    end if

                    if (size(state%vars(is)%values) /= 1) then
                        if (present(errmsg)) errmsg = 'State variable '//trim(vars(ivar)%name)//' must have exactly one value.'
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: Surface/bottom variable has wrong number of values for '//trim(vars(ivar)%name)
                        end if
                    end if

                    vars(ivar)%data_0d = state%vars(is)%values(1)

                case ('centre', 'interface')
                    if (.not. present(grid)) then
                        if (present(errmsg)) errmsg = 'Grid is required to initialise profile variable '//trim(vars(ivar)%name)
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: missing grid for profile variable '//trim(vars(ivar)%name)
                        end if
                    end if

                    if (.not. associated(vars(ivar)%data_1d)) then
                        if (present(errmsg)) errmsg = 'State variable '//trim(vars(ivar)%name)//' has no associated profile storage.'
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: unassociated profile storage for '//trim(vars(ivar)%name)
                        end if
                    end if

                    is_centre = trim(vars(ivar)%vert_coord) == 'centre'

                    if (is_centre) then
                        n_expected = grid%nz
                    else
                        n_expected = grid%nz + 1
                    end if

                    if (size(state%vars(is)%values) /= n_expected) then
                        if (present(errmsg)) errmsg = 'State variable '//trim(vars(ivar)%name)// &
                                                    ' has wrong profile size.'
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: wrong profile size for '//trim(vars(ivar)%name)
                        end if
                    end if

                    ! Validate depths
                    if (allocated(depth_work)) deallocate(depth_work)
                    if (allocated(value_work)) deallocate(value_work)
                    allocate(depth_work(n_expected), value_work(n_expected))
                    
                    if (allocated(tol)) deallocate(tol)
                    call default_profile_depth_tol(grid, tol, at_interfaces=.not. is_centre)

                    if (is_monotonic_decreasing(state%vars(is)%depth, mono_tol)) then
                        depth_work = state%vars(is)%depth
                        value_work = state%vars(is)%values

                    else if (is_monotonic_increasing(state%vars(is)%depth, mono_tol)) then
                        ! Flip arrays to go from bottom to surface
                        depth_work = state%vars(is)%depth(n_expected:1:-1)
                        value_work = state%vars(is)%values(n_expected:1:-1)

                    else
                        if (present(errmsg)) errmsg = 'Depth coordinate for state variable '//trim(vars(ivar)%name)// &
                                                      ' is not monotonic.'

                        if (allocated(depth_work)) deallocate(depth_work)
                        if (allocated(value_work)) deallocate(value_work)
                        if (allocated(tol)) deallocate(tol)
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: Depth coordinate needs to be monotonic, variable: '//trim(vars(ivar)%name)
                        end if
                    end if

                    ! Check that depths match
                    if (is_centre) then
                        depth_match = all(is_equal_array(depth_work, grid%z, abs_=tol))
                    else
                        depth_match = all(is_equal_array(depth_work, grid%z_w, abs_=tol))
                    end if
                    if (.not. depth_match) then
                        if (present(errmsg)) errmsg = 'Depth coordinate for state variable '//trim(vars(ivar)%name)// &
                                                      ' does not match model grid.'
                        
                        if (allocated(depth_work)) deallocate(depth_work)
                        if (allocated(value_work)) deallocate(value_work)
                        if (allocated(tol)) deallocate(tol)
                        if (present(ok)) then
                            ok = .false.
                            return
                        else
                            error stop 'set_initial_state: Depth does not match model grid for variable '//trim(vars(ivar)%name)
                        end if
                    end if
                    ! If everything is correct, then set the initial values
                    vars(ivar)%data_1d(:) = value_work(:)

                    if (allocated(depth_work)) deallocate(depth_work)
                    if (allocated(value_work)) deallocate(value_work)
                    if (allocated(tol)) deallocate(tol)
                case default
                    if (present(errmsg)) errmsg = 'Unsupported vertical coordinate for state variable '//trim(vars(ivar)%name)// &
                                                  ': '//trim(vars(ivar)%vert_coord)
                    if (present(ok)) then
                        ok = .false.
                        return
                    else
                        error stop 'set_initial_state: unsupported vert_coord for '//trim(vars(ivar)%name)
                    end if
            end select

        end do
        if (present(ok)) ok = .true.
    end subroutine set_initial_state


    !========================================================================
    !               Internal
    !========================================================================
    subroutine validate_state_path(filename, ok, errmsg)
        character(*), intent(in)  :: filename
        logical,      intent(out) :: ok
        character(*), intent(out) :: errmsg

        logical :: exists, supported

        exists = .false.
        supported = .false.
        ok = .false.
        errmsg = ''

        if (len_trim(filename) == 0) then
            errmsg = 'Initial state file path is empty.'
            return
        end if

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            errmsg = 'Initial state file does not exist: '//trim(filename)
            return
        end if

        supported = path_has_extension(filename, 'csv') .or. &
                    path_has_extension(filename, 'txt') .or. &
                    path_has_extension(filename, 'dat')


        if (.not. supported) then
            errmsg = 'Unsupported initial state file extension: '//trim(filename)// &
                     '. Expected .csv, .txt, or .dat. NetCDF restart files are not supported yet.'
            return
        end if

        ok = .true.
    end subroutine validate_state_path

    ! This is only for tidy formatted text files
    subroutine validate_state_columns(table, ok, errmsg)
        type(TextTable), intent(in)  :: table
        logical,         intent(out) :: ok
        character(*),    intent(out) :: errmsg

        integer :: i_name, i_depth, i_value

        ok = .false.
        errmsg = ''

        i_name  = table%find_column('name')
        i_depth = table%find_column('depth')
        i_value = table%find_column('value')

        if (i_name <= 0) then
            errmsg = 'Initial state file is missing required column: name'
            return
        end if

        if (i_depth <= 0) then
            errmsg = 'Initial state file is missing required column: depth'
            return
        end if

        if (i_value <= 0) then
            errmsg = 'Initial state file is missing required column: value'
            return
        end if

        ok = .true.
    end subroutine validate_state_columns

    subroutine build_state_from_tidy_table(table, state, ok, errmsg)
        type(TextTable), intent(in)    :: table
        type(StateData), intent(inout) :: state
        logical,         intent(out)   :: ok
        character(*),    intent(out)   :: errmsg

        integer :: i_name, i_depth, i_value, i_vert
        integer :: irow, ivar, nvars
        integer, allocatable :: counts(:), pos(:)
        character(len=NAME_LEN), allocatable :: names(:)
        character(len=16),       allocatable :: vert_coords(:)
        logical,                 allocatable :: has_vert_coords(:)
        character(len=NAME_LEN) :: vname
        character(len=16)       :: vcoord
        real(rk) :: z, val
        integer :: ios

        ok = .false.
        errmsg = ''

        i_name  = table%find_column('name')
        i_depth = table%find_column('depth')
        i_value = table%find_column('value')
        i_vert = table%find_column('vert_coord')

        nvars = 0
        allocate(names(table%nrow))
        allocate(counts(table%nrow))
        allocate(vert_coords(table%nrow))
        allocate(has_vert_coords(table%nrow))

        names  = ''
        counts = 0
        vert_coords = ''
        has_vert_coords = .false.

        ! First pass: identify variables and count rows.
        do irow = 1, table%nrow
            vname = trim(adjustl(table%values(irow, i_name)))
            vcoord = ''
            if (i_vert > 0) vcoord = trim(adjustl(table%values(irow, i_vert)))

            if (len_trim(vname) == 0) then
                errmsg = 'Empty variable name in initial state file at line '// &
                        trim(inttostr(table%line_numbers(irow)))//'.'
                return
            end if

            ivar = find_state_name(names, nvars, vname)

            if (ivar == 0) then
                nvars = nvars + 1
                names(nvars) = trim(vname)
                counts(nvars) = 1

                vert_coords(nvars) = trim(vcoord)
                has_vert_coords(nvars) = (i_vert > 0)
            else
                if (has_vert_coords(ivar)) then
                    if (trim(vert_coords(ivar)) /= trim(vcoord)) then
                        errmsg = 'Variable '//trim(vname)// &
                                ' appears with multiple vert_coord values.'
                        return
                    end if
                end if
                counts(ivar) = counts(ivar) + 1
            end if

        end do

        if (nvars == 0) then
            errmsg = 'Initial state file contains no data rows.'
            return
        end if

        allocate(state%vars(nvars))
        state%nvars = nvars

        do ivar = 1, nvars
            state%vars(ivar)%name = trim(names(ivar))
            state%vars(ivar)%vert_coord = trim(vert_coords(ivar))
            state%vars(ivar)%has_vert_coord = has_vert_coords(ivar)
            state%vars(ivar)%nvals = counts(ivar)

            allocate(state%vars(ivar)%depth(counts(ivar)))
            allocate(state%vars(ivar)%values(counts(ivar)))
        end do

        allocate(pos(nvars))
        pos = 0

        ! Second pass: fill variable-wise arrays.
        do irow = 1, table%nrow
            vname = trim(adjustl(table%values(irow, i_name)))
            ivar = find_state_name(names, nvars, vname)

            read(table%values(irow, i_depth), *, iostat=ios) z
            if (ios /= 0) then
                errmsg = 'Invalid depth value for variable '//trim(vname)// &
                        ' at line '//trim(inttostr(table%line_numbers(irow)))// &
                        ': '//trim(table%values(irow, i_depth))
                return
            end if
            if (state%vars(ivar)%nvals > 1 .and. z < 0.0_rk) then    
                errmsg = 'Negative depth for profile variable '//trim(vname)// &
                        ' at line '//trim(inttostr(table%line_numbers(irow)))//'.'
                return
            end if

            read(table%values(irow, i_value), *, iostat=ios) val
            if (ios /= 0) then
                errmsg = 'Invalid state value for variable '//trim(vname)// &
                        ' at line '//trim(inttostr(table%line_numbers(irow)))// &
                        ': '//trim(table%values(irow, i_value))
                return
            end if

            pos(ivar) = pos(ivar) + 1
            state%vars(ivar)%depth(pos(ivar))  = z
            state%vars(ivar)%values(pos(ivar)) = val
        end do

        ok = .true.
        if (allocated(names)) deallocate(names)
        if (allocated(counts)) deallocate(counts)
        if (allocated(vert_coords)) deallocate(vert_coords)
        if (allocated(has_vert_coords)) deallocate(has_vert_coords)
        if (allocated(pos)) deallocate(pos)
    end subroutine build_state_from_tidy_table

    integer function find_state_name(names, nnames, name) result(idx)
        character(len=*), intent(in) :: names(:)
        integer,          intent(in) :: nnames
        character(*),     intent(in) :: name

        integer :: i

        idx = 0

        do i = 1, nnames
            if (to_upper(trim(adjustl(names(i)))) == to_upper(trim(adjustl(name)))) then
                idx = i
                return
            end if
        end do
    end function find_state_name

end module state_loader