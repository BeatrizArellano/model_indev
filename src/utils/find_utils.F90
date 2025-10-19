! Useful simple functions to find elements in arrays
module find_utils
  use precision_types, only: rk    ! Importing real64 and int64

  implicit none
  private

  public :: find_name, has_name
  public :: argmin_abs_vec

contains

    ! Returns the index of the name to find in a vector of strings
    integer function find_name(list, name) result(idx)
        character(len=:), allocatable, intent(in) :: list(:)
        character(*), intent(in) :: name
        integer :: j
        idx = -1
        do j = 1, size(list)
            if (trim(adjustl(list(j))) == trim(adjustl(name))) then
                idx = j
                return
            end if
        end do
    end function find_name

    ! True if `name` matches any trimmed element of list (case-sensitive).
    !  list is OPTIONAL and ALLOCATABLE
    !  Use when the arrays may be unallocated
    logical function has_name(list, name)
        character(*), dimension(:), intent(in), optional :: list  ! allow optional
        character(*),               intent(in) :: name
        integer :: i
        has_name = .false.
        if (.not. present(list)) return
        do i=1, size(list)
            if (trim(adjustl(list(i))) == trim(adjustl(name))) then
                has_name = .true.; return
            end if
        end do
    end function has_name

    integer function argmin_abs_vec(v) result(idx)
      real(rk), intent(in) :: v(:)
      real(rk) :: best, cur
      integer :: i
      idx = 1; best = abs(v(1))
      do i=2, size(v)
         cur = abs(v(i))
         if (cur < best) then
            best = cur
            idx = i
         end if
      end do
   end function argmin_abs_vec






end module find_utils