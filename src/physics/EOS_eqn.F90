! eos.f90
module EOS_eqns
  use precision_types, only: rk
  implicit none
  private
  public :: eos_density, eos_mode_type

  enum, bind(c)
    enumerator :: EOS_UNESCO_1983 = 1
    enumerator :: EOS_LINEAR      = 2
    enumerator :: EOS_TEOS10      = 3   ! placeholder (via GSW later)
  end enum
  integer :: eos_mode_type = EOS_UNESCO_1983

contains

  pure elemental function eos_density(T, S) result(rho)  ! T[°C], S[PSU], p[dbar], rho[kg m-3]
    real(rk), intent(in)  :: T, S
    real(rk) :: rho
    ! Parameters for EOS_LINEAR
    real(rk), parameter :: rho0=1027._rk, T0=10._rk, S0=35._rk
    real(rk), parameter :: alpha=0.2_rk, beta=0.8_rk   ! kg m-3 per °C / PSU (rough)
    ! Parameters for EOS_UNESCO_1983
    real(rk) :: a0,a1,a2,a3,a4,a5, b0,b1,b2, c0,c1,c2
    select case (eos_mode_type)
        case (EOS_LINEAR)
            ! Simple linearized EOS about (T0,S0) for quick testing
            rho = rho0 - alpha*(T-T0) + beta*(S-S0)
        case (EOS_UNESCO_1983)
            ! UNESCO 1983 polynomial in T,S at atmospheric pressure (p ignored)
            ! (compact form)            
            a0=999.842594_rk; a1=6.793952e-2_rk; a2=-9.095290e-3_rk
            a3=1.001685e-4_rk; a4=-1.120083e-6_rk; a5=6.536332e-9_rk
            b0=8.24493e-1_rk;  b1=-4.0899e-3_rk;   b2=7.6438e-5_rk
            c0=-5.72466e-3_rk; c1=1.0227e-4_rk;    c2=-1.6546e-6_rk
            rho = (a0 + a1*T + a2*T*T + a3*T**3 + a4*T**4 + a5*T**5) &
                + (b0 + b1*T + b2*T*T)*S + (c0 + c1*T + c2*T*T)*S*sqrt(S)
            ! To add pressure correction later
    end select
  end function eos_density



end module EOS_eqns
