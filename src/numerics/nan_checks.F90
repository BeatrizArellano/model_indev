! src/utils/nan_checks.F90
module nan_checks
  use iso_fortran_env, only: error_unit
  use precision_types, only: rk
  use physics_types,   only: PhysicsState
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private

  public :: check_nan_physics

contains

    !-------------------------------------------------------------------
    ! Check for NaNs in physics-relevant fields.
    !
    ! 'where' is a label to know at what stage it was called, e.g.
    !   'after TURBULENCE_ke' or 'end of solve_physics step'.
    ! If a NaN is found, it prints a message to stderr and stops
    ! with a descriptive message.
    !-------------------------------------------------------------------
    subroutine check_nan_physics(PS, where)
      type(PhysicsState), intent(in) :: PS
      character(len=*), intent(in), optional :: where

      character(len=128) :: loc

      if (present(where)) then
        loc = trim(where)
      else
        loc = 'physics'
      end if

      call check_field_1d('temp',   Ps%temp,   loc)
      call check_field_1d('sal',    Ps%sal,    loc)
      call check_field_1d('rho',    PS%rho,    loc)
      call check_field_1d('velx',   PS%velx,   loc)
      call check_field_1d('vely',   PS%vely,   loc)
      call check_field_1d('Kz',     PS%Kz,     loc)
      call check_field_1d('Nz',     PS%Nz,     loc)
      call check_field_1d('tke',    PS%tke,    loc)
      call check_field_1d('eps',    PS%eps,    loc)
      call check_field_1d('Lscale', PS%Lscale, loc)
      call check_field_1d('NN',     PS%NN,     loc)
      call check_field_1d('SS',     PS%SS,     loc)
      call check_field_1d('Ri',     PS%Ri,     loc)

      ! Scalars
      call check_scalar('tau_x',  PS%tau_x,  loc)
      call check_scalar('tau_y',  PS%tau_y,  loc)
      call check_scalar('u_taus', PS%u_taus, loc)
      call check_scalar('u_taub', PS%u_taub, loc)
      call check_scalar('z0s',    PS%z0s,    loc)
      call check_scalar('z0b',    PS%z0b,    loc)
    end subroutine check_nan_physics


    !-------------------------------------------------------------------
    !Check a rank-1 real(rk) field for NaNs.
    ! Works for any lower/upper bounds (0..N or 1..N, etc.).
    !-------------------------------------------------------------------
    subroutine check_field_1d(name, a, loc)
      character(len=*), intent(in) :: name
      real(rk),         intent(in) :: a(:)
      character(len=*), intent(in) :: loc

      integer :: i, lb, ub

      lb = lbound(a, 1)
      ub = ubound(a, 1)

      do i = lb, ub
        if (ieee_is_nan(a(i))) then
          write(error_unit,'(A,": NaN in ",A," at index ",I0)') trim(loc), trim(name), i
          stop 'NaN detected in physics state'
        end if
      end do
    end subroutine check_field_1d

    subroutine check_scalar(name, x, loc)
        character(len=*), intent(in) :: name
        real(rk),         intent(in) :: x
        character(len=*), intent(in) :: loc

        if (ieee_is_nan(x)) then
            write(error_unit,'(A,": NaN in ",A)') trim(loc), trim(name)
            stop 'NaN detected in scalar physics variable'
        end if
    end subroutine check_scalar


end module nan_checks
