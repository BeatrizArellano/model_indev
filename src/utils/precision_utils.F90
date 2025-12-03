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
  pure elemental logical function is_equal(a, b, rel, abs_) result(eq)
    !! Robust floating-point comparison: true if |a-b| <= max(abs_tol, rel_tol*scale)
    real(rk), intent(in)           :: a, b
    real(rk), intent(in), optional :: rel, abs_
    real(rk) :: rtol, atol, scale

    ! Avoid MERGE with OPTIONAL – some compilers miscompile that pattern
    if (present(rel)) then
       rtol = rel
    else
       rtol = rel_default
    end if

    if (present(abs_)) then
       atol = abs_
    else
       atol = abs_default
    end if
    scale = max(1.0_rk, abs(a), abs(b))
    eq = abs(a - b) <= max(atol, rtol*scale)
  end function is_equal

  pure elemental logical function is_unequal(a, b, rel, abs_) result(neq)
    real(rk), intent(in)           :: a, b
    real(rk), intent(in), optional :: rel, abs_
    neq = .not. is_equal(a, b, rel, abs_)
  end function is_unequal

  pure elemental logical function is_zero(a, rel, abs_) result(isz)
    real(rk), intent(in)           :: a
    real(rk), intent(in), optional :: rel, abs_
    isz = is_equal(a, 0.0_rk, rel, abs_)
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
