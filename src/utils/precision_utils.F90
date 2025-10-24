! file: src/utils/precision_utils.f90
module precision_utils
  use, intrinsic :: ieee_arithmetic
  use precision_types, only: rk
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
  private
  public :: is_equal, is_unequal, is_zero
  public :: rel_default, abs_default, round_to

  ! Defaults scale with machine epsilon for rk
  real(rk), parameter :: rel_default = 1.0e4_rk * epsilon(1.0_rk)  ! ≈2e-12 for rk
  real(rk), parameter :: abs_default = 1.0e2_rk * epsilon(1.0_rk)  ! ≈2e-14 for rk

  interface round_to
     module procedure round_to_r32
     module procedure round_to_r64
  end interface

contains

  ! Work on scalars and arrays (element-wise)
  pure elemental logical function is_equal(a, b, rel, abs_)
    !   Robust floating-point comparison: returns .true. if A and B are
    !   equal within a combined tolerance band (relative + absolute).
    !   Arguments:
    !     a, b   [in]  real(rk)  - values to compare.
    !     rel    [in]  real(rk), optional - relative tolerance (dimensionless).
    !                         If absent, uses rel_default.
    !     abs_   [in]  real(rk), optional - absolute tolerance (same units as a/b).
    !                         If absent, uses abs_default.
    real(rk), intent(in)           :: a, b
    real(rk), intent(in), optional :: rel, abs_
    real(rk) :: rtol, atol, scale
    rtol  = merge(rel,  rel_default, present(rel))
    atol  = merge(abs_, abs_default, present(abs_))
    scale = max(1.0_rk, abs(a), abs(b))
    is_equal = abs(a - b) <= max(atol, rtol*scale)
  end function is_equal

  pure elemental logical function is_unequal(a, b, rel, abs_)
  !   Negation of is_equal: returns .true. if A and B differ by more
  !     than the tolerance band (relative + absolute).
    real(rk), intent(in)           :: a, b
    real(rk), intent(in), optional :: rel, abs_
    is_unequal = .not. is_equal(a, b, rel, abs_)
  end function is_unequal

  pure elemental logical function is_zero(a, rel, abs_)
    ! Checks whether a value a is within tolerance of 0.0_rk.
    real(rk), intent(in)           :: a
    real(rk), intent(in), optional :: rel, abs_
    is_zero = is_equal(a, 0.0_rk, rel, abs_)
  end function is_zero


  pure elemental function round_to_r32(x, ndigits) result(r)
    real(real32), intent(in) :: x
    integer,      intent(in) :: ndigits
    real(real32)             :: r
    real(real32)             :: scale, ten, y

    ! Propagate NaN/Inf unchanged
    if (.not. ieee_is_finite(x)) then
       r = x
       return
    end if

    ten = real(10.0, kind=real32)

    if (ndigits == 0) then
       r = anint(x)                         ! ties-to-even under IEEE nearest (default)
    else if (ndigits > 0) then
       ! Round to +ndigits decimal places
       ! Guard extreme exponents: if scale over/underflows, just return x
       if (ndigits > 300) then
          r = x
       else
          scale = ten**ndigits
          y     = x * scale
          r     = anint(y) / scale
       end if
    else
       ! ndigits < 0: round to tens/hundreds/etc.
       if (ndigits < -300) then
          r = anint(x) * 0.0_real32    ! effectively 0 with sign of anint(x); avoids inf scale
       else
          scale = ten**(-ndigits)
          y     = x / scale
          r     = anint(y) * scale
       end if
    end if
  end function round_to_r32


  pure elemental function round_to_r64(x, ndigits) result(r)
      real(rk), intent(in) :: x
      integer,      intent(in) :: ndigits
      real(rk)             :: r
      real(rk)             :: scale, ten, y

      if (.not. ieee_is_finite(x)) then
        r = x
        return
      end if

      ten = real(10.0, kind=rk)

      if (ndigits == 0) then
        r = anint(x)
      else if (ndigits > 0) then
        if (ndigits > 300) then
            r = x
        else
            scale = ten**ndigits
            y     = x * scale
            r     = anint(y) / scale
        end if
      else
        if (ndigits < -300) then
            r = anint(x) * 0.0_rk
        else
            scale = ten**(-ndigits)
            y     = x / scale
            r     = anint(y) * scale
        end if
      end if
  end function round_to_r64

end module precision_utils
