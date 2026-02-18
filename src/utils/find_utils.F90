! Useful simple functions to find elements in arrays
module find_utils
  use precision_types, only: rk    ! Importing real64 and int64

  implicit none
  private

  public :: find_name, has_name
  public :: argmin_abs_vec, find_nearest_index

contains

    ! Returns the index of the element in array `a` closest to `x` (minimizes |a(i)−x|).
    ! Tie-break (equal distance): if prefer_min=true choose smaller a(i), else choose larger a(i).
    ! If still tied (equal values), keep the first occurrence.
    integer function find_nearest_index(a, x, prefer_min) result(idx)
        real(rk), intent(in) :: a(:)
        real(rk), intent(in) :: x
        logical, intent(in), optional :: prefer_min

        real(rk) :: best, cur
        logical  :: pmin
        integer  :: i

        pmin = .true.
        if (present(prefer_min)) pmin = prefer_min

        idx  = 1
        best = abs(a(1) - x)

        do i = 2, size(a)
            cur = abs(a(i) - x)

            if (cur < best) then
                best = cur
                idx  = i

            else if (cur == best) then
                if (pmin) then
                    if (a(i) < a(idx)) idx = i
                else
                    if (a(i) > a(idx)) idx = i
                end if
                ! If a(i) == a(idx), keep earlier index (stable)
            end if
        end do
    end function find_nearest_index



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