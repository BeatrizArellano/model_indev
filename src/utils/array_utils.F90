module array_utils
  use precision_types, only: rk
  implicit none
  private

  public :: reverse_array
  public :: str_to_real_vec

  ! Generic interface so reverse_array works for different types
  interface reverse_array
     module procedure reverse_real_1d
     module procedure reverse_int_1d
  end interface reverse_array

contains      

  subroutine str_to_real_vec(strings, values, ok, errmsg, label)
      character(*), intent(in)               :: strings(:)
      real(rk), allocatable, intent(out)     :: values(:)
      logical, intent(out), optional         :: ok
      character(*), intent(out), optional    :: errmsg
      character(*), intent(in), optional     :: label

      integer :: i
      integer :: ios

      logical :: lok
      character(len=512) :: lmsg

      lok  = .true.
      lmsg = ''

      allocate(values(size(strings)))

      do i = 1, size(strings)
          read(strings(i), *, iostat=ios) values(i)

          if (ios /= 0) then
              lok = .false.
              if (present(label)) then
                  write(lmsg,'(A,I0,A,A,A)') &
                      'Failed to parse value at index ', i, &
                      ' in column "', trim(label), '".'
              else
                  write(lmsg,'(A,I0,A)') &
                      'Failed to parse value at index ', i, '.'
              end if
              exit
          end if
      end do

      if (present(ok)) then
          ok = lok
      else
          if (.not. lok) error stop trim(lmsg)
      end if
      if (present(errmsg)) errmsg = trim(lmsg)

  end subroutine str_to_real_vec



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
