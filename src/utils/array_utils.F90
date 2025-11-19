module array_utils
  use precision_types, only: rk
  implicit none
  private

  public :: reverse_array

  ! Generic interface so reverse_array works for different types
  interface reverse_array
     module procedure reverse_real_1d
     module procedure reverse_int_1d
  end interface reverse_array

contains

  !-------------------------------------------------------------------
  ! Reverse a real(rk) 1D array.
  !
  !  src(:)  : input array (any bounds, not modified)
  !  dest(:) : output array (allocated here, same bounds as src)
  !
  pure subroutine reverse_real_1d(src, dest)
    real(rk), intent(in)              :: src(:)
    real(rk), allocatable, intent(out):: dest(:)
    integer :: lb, ub, i

    lb = lbound(src, 1)
    ub = ubound(src, 1)

    allocate(dest(lb:ub))

    do i = lb, ub
       dest(i) = src(ub - (i - lb))
    end do
  end subroutine reverse_real_1d

  !-------------------------------------------------------------------
  ! Reverse an integer 1D array.
  !
  !  src(:)  : input array (any bounds, not modified)
  !  dest(:) : output array (allocated here, same bounds as src)
  !
  pure subroutine reverse_int_1d(src, dest)
    integer, intent(in)               :: src(:)
    integer, allocatable, intent(out) :: dest(:)
    integer :: lb, ub, i

    lb = lbound(src, 1)
    ub = ubound(src, 1)

    allocate(dest(lb:ub))

    do i = lb, ub
       dest(i) = src(ub - (i - lb))
    end do
  end subroutine reverse_int_1d

end module array_utils
