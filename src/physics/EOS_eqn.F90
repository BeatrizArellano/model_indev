! eos.f90
module EOS_eqns
  use precision_types, only: rk
  implicit none
  private
  public :: eos_density, eos_mode_type

  enum, bind(c)
    enumerator :: EOS_UNESCO_1983 = 1
    enumerator :: EOS_LINEAR      = 2
    enumerator :: EOS_TEOS10      = 3   ! Not implemented yet
  end enum
  integer :: eos_mode_type = EOS_UNESCO_1983

contains

  pure elemental function eos_density(T, S) result(rho)
    !-----------------------------------------------------------------
    ! Compute seawater density rho [kg m-3] from temperature T [°C]
    ! and practical salinity S [PSU].
    !
    ! EOS_LINEAR      : simple linearised EOS about (T0, S0)
    ! EOS_UNESCO_1983 : UNESCO 1983 (EOS-80) density at p = 0 dbar
    !-----------------------------------------------------------------
    real(rk), intent(in) :: T, S
    real(rk)             :: rho

    ! Linear EOS parameters 
    real(rk), parameter :: rho0  = 1027._rk   ! reference density [kg m-3]
    real(rk), parameter :: T0    = 10._rk     ! reference temperature [°C]
    real(rk), parameter :: S0    = 35._rk     ! reference salinity [PSU]
    real(rk), parameter :: alpha = 0.2_rk     ! d(rho)/dT   [kg m-3 °C-1]
    real(rk), parameter :: beta  = 0.8_rk     ! d(rho)/dS   [kg m-3 PSU-1]

    ! UNESCO 1983 (EOS-80, p =0)
    real(rk) :: a0, a1, a2, a3, a4, a5
    real(rk) :: rho_w, A, B, C, sqrtS

    select case (eos_mode_type)
      case (EOS_LINEAR)
          ! Simple linearized EOS for quick testing
          rho = rho0 - alpha*(T - T0) + beta*(S - S0)

      case (EOS_UNESCO_1983)
          ! UNESCO 1983 (EOS-80) density at atmospheric pressure (p = 0 dbar)
          ! rho(S,T) = rho_w(T) + A(T)*S + B(T)*S^(3/2) + C*S^2

          ! Pure water density rho_w(T)
          a0 = 999.842594_rk
          a1 = 6.793952e-2_rk
          a2 = -9.095290e-3_rk
          a3 = 1.001685e-4_rk
          a4 = -1.120083e-6_rk
          a5 = 6.536332e-9_rk

          rho_w = a0 + a1*T + a2*T*T + a3*T**3 + a4*T**4 + a5*T**5

          ! A(T): coefficient for the linear-S term
          A =  8.24493e-1_rk          &
            - 4.0899e-3_rk*T         &
            + 7.6438e-5_rk*T*T       &
            - 8.2467e-7_rk*T**3      &
            + 5.3875e-9_rk*T**4

          ! B(T): coefficient for the S^(3/2) term
          B = -5.72466e-3_rk          &
            + 1.0227e-4_rk*T         &
            - 1.6546e-6_rk*T*T

          ! C: coefficient for the S^2 term (constant)
          C = 4.8314e-4_rk

          sqrtS = sqrt(S)
          rho = rho_w + A*S + B*S*sqrtS + C*S*S

      case default
          error stop 'eos_density: unknown eos_mode_type'

    end select

end function eos_density




end module EOS_eqns
