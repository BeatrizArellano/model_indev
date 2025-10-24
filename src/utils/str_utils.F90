! Useful functions for strings
module str_utils
  use precision_types, only: rk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  implicit none
  private

  public :: to_lower, replace_char_inplace
  public :: inttostr, realtostr
  public :: list_to_str, append_string
contains
 
    ! Convert string to lowercase
    ! Notes: Only 'A'..'Z' are mapped; other chars are unchanged.
    pure function to_lower(s) result(t)
        character(*), intent(in) :: s
        character(len=len(s)) :: t
        integer :: i, ich
        do i=1,len(s)
            ich = iachar(s(i:i))
            if (ich>=iachar('A') .and. ich<=iachar('Z')) then
                t(i:i) = achar(ich + 32)
            else
                t(i:i) = s(i:i)
            end if
        end do
    end function to_lower

    ! Replace all occurrences of a character in-place, returning new string
    ! Replace character a with b in string s
    subroutine replace_char_inplace(s, a, b)
        character(len=:), allocatable, intent(inout) :: s
        character(*), intent(in) :: a, b
        integer :: i
        do i=1, len(s)
            if (s(i:i) == a) s(i:i) = b
        end do
    end subroutine replace_char_inplace

    ! Converts an integer to string
    pure function inttostr(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: tmp
        write(tmp,'(I0)') i
        s = trim(tmp)
    end function inttostr

    ! Converts a real to string
    pure function realtostr(x, ndp) result(s)  
        real(rk), intent(in)           :: x
        integer,  intent(in), optional :: ndp
        character(len=:), allocatable  :: s

        character(len=64) :: tmp
        character(len=24) :: fmt
        integer :: p, n

        ! Special values first
        if (.not. ieee_is_finite(x)) then
            if (ieee_is_nan(x)) then
                s = 'NaN'
            else
                if (x > 0.0_rk) then
                s = 'Inf'
                else
                s = '-Inf'
                end if
            end if
            return
        end if

        if (present(ndp)) then
            ! Clamp to a sensible precision range
            p = max(0, min(15, ndp))      ! adjust 15 to your rk precision if needed
            write(fmt, '("(F0.", I0, ")")') p
            write(tmp, fmt) x
            s = adjustl(trim(tmp))

            ! Trim trailing zeros and a trailing '.' for a cleaner look
            if (index(s, '.') > 0) then
                n = len_trim(s)
                do while (n > 0 .and. s(n:n) == '0')
                n = n - 1
                end do
                if (n > 0 .and. s(n:n) == '.') n = n - 1
                s = s(:n)
            end if
        else
            ! Compact output with minimal width
            write(tmp, '(G0)') x
            s = adjustl(trim(tmp))
        end if
    end function realtostr

    ! Joins list elements in a string using sep as separator
    ! Example list_to_str([e1,..en],',') -> 'e1,e2,...,en'
    pure function list_to_str(arr, sep) result(s)
        character(len=:), allocatable, intent(in) :: arr(:)
        character(*), intent(in) :: sep
        character(len=:), allocatable :: s
        integer :: i
        if (.not. allocated(arr) .or. size(arr) == 0) then
            s = ''; return
        end if
        s = trim(arr(1))
        do i = 2, size(arr)
            s = s // sep // trim(arr(i))
        end do
    end function list_to_str

    ! Appends trimmed item to allocatable character array,
    ! allocating/expanding as needed while preserving existing contents.
    subroutine append_string(arr, item)
        character(:), allocatable, intent(inout) :: arr(:)
        character(*),               intent(in)    :: item
        integer :: n
        if (.not. allocated(arr)) then
            allocate(character(len=len_trim(item)) :: arr(1))
            arr(1) = trim(item)
        else
            n = size(arr)
            call extend_char_array(arr, n+1, len_trim(item))
            arr(n+1) = trim(item)
        end if
    end subroutine append_string

    !----------------------------------
    ! Local Helpers
    !------------------------------------

    ! Resize allocatable character array to length new_n and element LEN=new_len, preserving existing elements.
    subroutine extend_char_array(arr, new_n, new_len)
        character(:), allocatable, intent(inout) :: arr(:)
        integer,                    intent(in)    :: new_n, new_len
        character(:), allocatable :: tmp(:)
        integer :: i, old_n, L
        if (.not. allocated(arr)) then
            allocate(character(len=new_len) :: arr(new_n))
            return
        end if
        old_n = size(arr); L = len(arr(1))
        allocate(character(len=max(L,new_len)) :: tmp(new_n))
        tmp = ''
        do i=1, old_n
            tmp(i)(1:min(len(arr(i)),len(tmp(i)))) = arr(i)(1:min(len(arr(i)),len(tmp(i))))
        end do
        call move_alloc(tmp, arr)
    end subroutine extend_char_array

end module str_utils