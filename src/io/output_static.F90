module output_static
    use precision_types, only: rk
    implicit none
    private 
    public :: StaticProfile
    public :: append_static_profile
    public :: init_static_profile
    public :: clear_static_profiles 
    !===========================================================
    ! Type definition
    !===========================================================
    type :: StaticProfile
       character(len=:), allocatable :: name
       character(len=:), allocatable :: long_name
       character(len=:), allocatable :: units
       character(len=:), allocatable :: standard_name   
       character(len=16) :: vert_coord = 'centre'  ! 'centre' or 'interface'    
       logical :: has_min     = .false.
       logical :: has_max     = .false.
       logical :: has_missing = .false. 
       real(rk) :: valid_min     = 0.0_rk
       real(rk) :: valid_max     = 0.0_rk
       real(rk) :: missing_value = 0.0_rk   
       real(rk), allocatable :: data_1d(:)
    end type StaticProfile  
contains

    !===========================================================
    ! Initialise a profile 
    !===========================================================
    subroutine init_static_profile(p, name, long_name, units, profile_data, &
                                   standard_name, vert_coord)
        type(StaticProfile), intent(out) :: p
        character(len=*),    intent(in)  :: name
        character(len=*),    intent(in)  :: long_name
        character(len=*),    intent(in)  :: units
        real(rk),            intent(in)  :: profile_data(:)
        character(len=*),    intent(in), optional :: standard_name
        character(len=*),    intent(in), optional :: vert_coord

        p%has_min       = .false.
        p%has_max       = .false.
        p%has_missing   = .false.
        p%valid_min     = 0.0_rk
        p%valid_max     = 0.0_rk
        p%missing_value = 0.0_rk
        p%standard_name = ''   
        
        if (len_trim(name) == 0) then
            error stop 'init_static_profile: profile name cannot be empty.'
        end if

        p%name      = trim(name)
        p%long_name = trim(long_name)
        p%units     = trim(units)

        if (present(standard_name)) then
            if (len_trim(standard_name) > 0) p%standard_name = trim(standard_name)
        end if
        if (present(vert_coord)) then
            select case (trim(adjustl(vert_coord)))
            case ('centre', 'interface')
                p%vert_coord = trim(adjustl(vert_coord))
            case default
                error stop 'init_static_profile: vert_coord must be "centre" or "interface".'
            end select
        else
            p%vert_coord = 'centre'
        end if

        allocate(p%data_1d(size(profile_data)))
        p%data_1d = profile_data
    end subroutine init_static_profile


    !===========================================================
    ! Append a profile to a list
    !===========================================================
    subroutine append_static_profile(list, prof, allow_replace)
        type(StaticProfile), allocatable, intent(inout) :: list(:)
        type(StaticProfile),              intent(in)    :: prof
        logical, optional,                intent(in)    :: allow_replace

        type(StaticProfile), allocatable :: tmp(:)
        integer :: n, i
        logical :: replace

        replace = .false.
        if (present(allow_replace)) replace = allow_replace

        !---------------------------------------------------------
        ! If list is not allocated then initialise
        !---------------------------------------------------------
        if (.not. allocated(list)) then
            allocate(list(1))
            call copy_profile(list(1), prof)
            return
        end if

        !---------------------------------------------------------
        ! Check for duplicate names
        !---------------------------------------------------------
        do i = 1, size(list)
            if (trim(adjustl(list(i)%name)) == trim(adjustl(prof%name))) then
                if (replace) then
                    call copy_profile(list(i), prof)
                    return
                else
                    error stop 'append_static_profile: duplicate profile name: '//trim(prof%name)
                end if
            end if
        end do

        !---------------------------------------------------------
        ! Append (reallocate)
        !---------------------------------------------------------
        n = size(list)
        allocate(tmp(n+1))

        ! Copy existing
        do i = 1, n
            call copy_profile(tmp(i), list(i))
        end do

        ! Add new
        call copy_profile(tmp(n+1), prof)

        call move_alloc(tmp, list)

    end subroutine append_static_profile


    !===========================================================
    ! Clear list
    !===========================================================
    subroutine clear_static_profiles(list)
        type(StaticProfile), allocatable, intent(inout) :: list(:)
        if (allocated(list)) deallocate(list)
    end subroutine clear_static_profiles


    !===========================================================
    ! Internal: deep copy
    !===========================================================
    subroutine copy_profile(dest, src)
        type(StaticProfile), intent(out) :: dest
        type(StaticProfile), intent(in)  :: src

        dest%name          = src%name
        dest%long_name     = src%long_name
        dest%units         = src%units
        if (allocated(src%standard_name)) dest%standard_name = src%standard_name
        dest%vert_coord    = src%vert_coord

        dest%has_min       = src%has_min
        dest%has_max       = src%has_max
        dest%has_missing   = src%has_missing

        dest%valid_min     = src%valid_min
        dest%valid_max     = src%valid_max
        dest%missing_value = src%missing_value

        if (allocated(src%data_1d)) then
            allocate(dest%data_1d(size(src%data_1d)))
            dest%data_1d = src%data_1d
        end if

    end subroutine copy_profile

end module output_static