module variable_registry
    use precision_types, only: rk
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    private

    type :: VarMetadata
        character(len=:), allocatable :: name
        character(len=:), allocatable :: standard_name   ! CF standard_name
        character(len=:), allocatable :: long_name
        character(len=:), allocatable :: units

        character(len=16) :: vert_coord                  ! none, centre, interface, surface, bottom
        integer           :: n_space_dims                ! 0 = scalar, 1 = profile (depth)       

        logical           :: output = .false.            ! Whether to include this variable in the output file

         ! Optional value metadata
        logical           :: has_min      = .false.
        logical           :: has_max      = .false.
        logical           :: has_missing  = .false.
        real(rk)          :: valid_min    = 0._rk
        real(rk)          :: valid_max    = 0._rk
        real(rk)          :: missing_value = 0._rk

        real(rk), pointer :: data_1d(:) => null()
        real(rk), pointer :: data_0d    => null()
    end type VarMetadata

    public :: VarMetadata
    public :: register_variable
    public :: output_all_variables, enable_named_variables

    interface register_variable
        module procedure register_variable_0d
        module procedure register_variable_1d
    end interface register_variable

contains
    subroutine register_variable_1d(registry, name, long_name, units, vert_coord, n_space_dims, data_1d, &
                                    standard_name, minimum, maximum, missing_value)
        type(VarMetadata), allocatable, intent(inout) :: registry(:)
        character(len=*),                intent(in)    :: name, long_name, units, vert_coord
        character(len=*),      optional, intent(in)    :: standard_name
        integer,                         intent(in)    :: n_space_dims
        real(rk),         target,        intent(in)    :: data_1d(:)
        real(rk),         optional,      intent(in)    :: minimum, maximum, missing_value


        type(VarMetadata), allocatable :: tmp(:)
        integer :: n

        if (.not. allocated(registry)) then
            allocate(registry(1))
            n = 1
        else
            n = size(registry) + 1
            allocate(tmp(n))
            tmp(1:n-1) = registry
            call move_alloc(tmp, registry)
        end if

        registry(n)%name      = trim(name)

        if (present(standard_name)) then
            registry(n)%standard_name = trim(standard_name)
        else
            registry(n)%standard_name = ''   ! no CF name
        end if

        registry(n)%long_name = trim(long_name)
        registry(n)%units     = trim(units)
        registry(n)%vert_coord  = trim(vert_coord)
        registry(n)%n_space_dims = n_space_dims
        registry(n)%output   = .false.
        ! Value metadata
        if (present(minimum) .and. .not. ieee_is_nan(minimum)) then
            registry(n)%has_min   = .true.
            registry(n)%valid_min = minimum
        else
            registry(n)%has_min   = .false.
        end if

        if (present(maximum) .and. .not. ieee_is_nan(maximum)) then
            registry(n)%has_max   = .true.
            registry(n)%valid_max = maximum
        else
            registry(n)%has_max   = .false.
        end if

        if (present(missing_value).and. .not. ieee_is_nan(missing_value)) then
            registry(n)%has_missing   = .true.
            registry(n)%missing_value = missing_value
        else
            registry(n)%has_missing   = .false.
        end if

        nullify(registry(n)%data_0d)
        nullify(registry(n)%data_1d)
        ! Pointers to data
        registry(n)%data_1d  => data_1d
    end subroutine register_variable_1d

    subroutine register_variable_0d(registry, name, long_name, units, vert_coord, n_space_dims, data_0d, &
                                    standard_name, minimum, maximum, missing_value)
        type(VarMetadata), allocatable, intent(inout) :: registry(:)
        character(len=*),                intent(in)    :: name, long_name, units, vert_coord
        character(len=*),      optional, intent(in)    :: standard_name
        integer,                         intent(in)    :: n_space_dims
        real(rk),         target,        intent(in)    :: data_0d
        real(rk),         optional,      intent(in)    :: minimum, maximum, missing_value

        type(VarMetadata), allocatable :: tmp(:)
        integer :: n

        if (.not. allocated(registry)) then
            allocate(registry(1))
            n = 1
        else
            n = size(registry) + 1
            allocate(tmp(n))
            tmp(1:n-1) = registry
            call move_alloc(tmp, registry)
        end if

        registry(n)%name      = trim(name)
        if (present(standard_name)) then
            registry(n)%standard_name = trim(standard_name)
        else
            registry(n)%standard_name = ''
        end if
        registry(n)%long_name = trim(long_name)
        registry(n)%units     = trim(units)
        registry(n)%vert_coord   = trim(vert_coord)
        registry(n)%n_space_dims = n_space_dims
        registry(n)%output   = .false.

        ! Value metadata
        if (present(minimum) .and. .not. ieee_is_nan(minimum)) then
            registry(n)%has_min   = .true.
            registry(n)%valid_min = minimum
        else
            registry(n)%has_min   = .false.
        end if

        if (present(maximum) .and. .not. ieee_is_nan(maximum)) then
            registry(n)%has_max   = .true.
            registry(n)%valid_max = maximum
        else
            registry(n)%has_max   = .false.
        end if

        if (present(missing_value).and. .not. ieee_is_nan(missing_value)) then
            registry(n)%has_missing   = .true.
            registry(n)%missing_value = missing_value
        else
            registry(n)%has_missing   = .false.
        end if

        nullify(registry(n)%data_0d)
        registry(n)%data_0d  => data_0d
        nullify(registry(n)%data_1d)
    end subroutine register_variable_0d

    subroutine output_all_variables(registry)
        type(VarMetadata), allocatable, intent(inout) :: registry(:)
        integer :: i
        if (.not. allocated(registry)) return
        do i = 1, size(registry)
            registry(i)%output = .true.
        end do
    end subroutine output_all_variables

    subroutine enable_named_variables(registry, list, n_list)
        type(VarMetadata), allocatable, intent(inout) :: registry(:)
        character(len=*),           intent(in)    :: list(:)
        integer,                    intent(in)    :: n_list
        integer :: i, j

        if (.not. allocated(registry)) return

        do i = 1, size(registry)
            registry(i)%output = .false.
        end do

        do j = 1, n_list
            do i = 1, size(registry)
                if (trim(registry(i)%name) == trim(list(j))) registry(i)%output = .true.
            end do
        end do
    end subroutine enable_named_variables

end module variable_registry
