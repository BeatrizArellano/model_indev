module tridiagonal
  use precision_types, only: rk
  implicit none
  private
  public :: TridiagCoeff, init_tridiag, reset_tridiag, solve_tridiag

  type :: TridiagCoeff
     real(rk), allocatable :: au(:), bu(:), cu(:), du(:)  ! coefficients
     real(rk), allocatable :: ru(:), qu(:)                ! scratch
  end type TridiagCoeff

contains

  subroutine init_tridiag(coeff, N)
    type(TridiagCoeff), intent(inout) :: coeff
    integer,            intent(in)    :: N
    if (.not. allocated(coeff%au)) allocate(coeff%au(0:N))
    if (.not. allocated(coeff%bu)) allocate(coeff%bu(0:N))
    if (.not. allocated(coeff%cu)) allocate(coeff%cu(0:N))
    if (.not. allocated(coeff%du)) allocate(coeff%du(0:N))
    if (.not. allocated(coeff%ru)) allocate(coeff%ru(0:N))
    if (.not. allocated(coeff%qu)) allocate(coeff%qu(0:N))
    call reset_tridiag(coeff)
  end subroutine init_tridiag

  pure subroutine reset_tridiag(coeff)
    type(TridiagCoeff), intent(inout) :: coeff
    coeff%au = 0.0_rk; coeff%bu = 0.0_rk; coeff%cu = 0.0_rk; coeff%du = 0.0_rk
    coeff%ru = 0.0_rk; coeff%qu = 0.0_rk
  end subroutine reset_tridiag

  ! Thomas algorithm to solve a tridiagonal matrix
  subroutine solve_tridiag(N, fi, lt, coeff, value)
    integer,            intent(in)  :: N, fi, lt
    type(TridiagCoeff), intent(inout) :: coeff
    real(rk),           intent(out) :: value(:)

    integer :: i
    real(rk) :: denom
    real(rk), parameter :: eps = 1.0e-30_rk

    ! Forward sweep
    !denom = merge(coeff%bu(lt), sign(1.0_rk, coeff%bu(lt))*eps, abs(coeff%bu(lt))>eps) ! Ensure the diagonal s big enough
    coeff%ru(lt) = coeff%au(lt) / coeff%bu(lt)
    coeff%qu(lt) = coeff%du(lt) / coeff%bu(lt)

    do i = lt-1, fi+1, -1
       denom        = coeff%bu(i) - coeff%cu(i)*coeff%ru(i+1)
       if (abs(denom) <= eps) denom = sign(eps, denom)  ! safeguard
       coeff%ru(i) = coeff%au(i) / denom
       coeff%qu(i) = (coeff%du(i) - coeff%cu(i)*coeff%qu(i+1)) / denom
    end do

    denom      = coeff%bu(fi) - coeff%cu(fi)*coeff%ru(fi+1)
    if (abs(denom) <= eps) denom = sign(eps, denom)
    coeff%qu(fi) = (coeff%du(fi) - coeff%cu(fi)*coeff%qu(fi+1)) / denom

    ! Back substitution
    value(fi) = coeff%qu(fi)
    do i = fi+1, lt
       value(i) = coeff%qu(i) - coeff%ru(i)*value(i-1)
    end do
  end subroutine solve_tridiag

end module tridiagonal
