module tridiagonal
  use precision_types, only: rk
  implicit none
  private
  public :: TridiagCoeff, init_tridiag, reset_tridiag, solve_tridiag, clear_tridiag

  type :: TridiagCoeff
     real(rk), allocatable :: au(:), bu(:), cu(:), du(:)  ! coefficients
     real(rk), allocatable :: ru(:), qu(:)                ! scratch
  end type TridiagCoeff

contains

  subroutine init_tridiag(coeff, N)
    type(TridiagCoeff), intent(inout) :: coeff
    integer,            intent(in)    :: N
    if (.not. allocated(coeff%au)) allocate(coeff%au(N))
    if (.not. allocated(coeff%bu)) allocate(coeff%bu(N))
    if (.not. allocated(coeff%cu)) allocate(coeff%cu(N))
    if (.not. allocated(coeff%du)) allocate(coeff%du(N))
    if (.not. allocated(coeff%ru)) allocate(coeff%ru(N))
    if (.not. allocated(coeff%qu)) allocate(coeff%qu(N))
    call reset_tridiag(coeff)
  end subroutine init_tridiag

  pure subroutine reset_tridiag(coeff)
    type(TridiagCoeff), intent(inout) :: coeff
    coeff%au = 0.0_rk; coeff%bu = 0.0_rk; coeff%cu = 0.0_rk; coeff%du = 0.0_rk
    coeff%ru = 0.0_rk; coeff%qu = 0.0_rk
  end subroutine reset_tridiag

  ! Thomas algorithm to solve a tridiagonal matrix
  subroutine solve_tridiag(fi, lt, coeff, values)
      integer,            intent(in)  :: fi, lt
      type(TridiagCoeff), intent(inout) :: coeff
      real(rk),           intent(out) :: values(:)

      integer :: i
      real(rk) :: denom
      real(rk), parameter :: eps = 1.0e-30_rk

      ! Forward sweep
      coeff%ru(lt) = coeff%au(lt) / coeff%bu(lt)
      coeff%qu(lt) = coeff%du(lt) / coeff%bu(lt)

      do i = lt-1, fi+1, -1
        denom        = coeff%bu(i) - coeff%cu(i)*coeff%ru(i+1)
        !if (abs(denom) <= eps) denom = sign(eps, denom)  ! safeguard
        coeff%ru(i) = coeff%au(i) / denom
        coeff%qu(i) = (coeff%du(i) - coeff%cu(i)*coeff%qu(i+1)) / denom
      end do

      denom      = coeff%bu(fi) - coeff%cu(fi)*coeff%ru(fi+1)
      if (abs(denom) <= eps) denom = sign(eps, denom)
      coeff%qu(fi) = (coeff%du(fi) - coeff%cu(fi)*coeff%qu(fi+1)) / denom

      ! Back substitution
      values(fi) = coeff%qu(fi)
      do i = fi+1, lt
        values(i) = coeff%qu(i) - coeff%ru(i)*values(i-1)
      end do

  end subroutine solve_tridiag

  ! Clears space in memory by deallocating arrays
  subroutine clear_tridiag(coeff)
      type(TridiagCoeff), intent(inout) :: coeff

      if (allocated(coeff%au)) deallocate(coeff%au)
      if (allocated(coeff%bu)) deallocate(coeff%bu)
      if (allocated(coeff%cu)) deallocate(coeff%cu)
      if (allocated(coeff%du)) deallocate(coeff%du)
      if (allocated(coeff%ru)) deallocate(coeff%ru)
      if (allocated(coeff%qu)) deallocate(coeff%qu)
  end subroutine clear_tridiag

end module tridiagonal
