module pressure
  use precision_types, only: rk
  use physics_params,  only: gravity, p_atm_dbar
  implicit none
  private
  public :: compute_pressure



contains

  !---------------------------------------------------------------
  ! Compute hydrostatic pressure at layer centres [dbar], including
  ! atmospheric pressure.
  !
  ! Inputs:
  !   dz(1:N)    : layer thicknesses [m], positive downward
  !   rho(1:N)   : density at layer centres [kg m-3]
  !   p_air_dbar : optional atmospheric pressure [dbar]
  !
  ! Output:
  !   pres(1:N)  : absolute pressure at centres [dbar]
  !
  ! Convention:
  !   - index 1 = bottom layer, index N = surface layer
  !---------------------------------------------------------------
  pure subroutine compute_pressure(dz, rho, pres, p_air_dbar)
    real(rk),           intent(in)  :: dz(:), rho(:)
    real(rk),           intent(out) :: pres(:)
    real(rk), optional, intent(in)  :: p_air_dbar

    integer  :: N, i
    real(rk) :: p0

    N = size(rho)

    ! Atmospheric pressure at the surface [dbar]
    p0 = p_atm_dbar
    if (present(p_air_dbar)) p0 = p_air_dbar

    ! Build integral of rho*dz in "pres" (kg m-2), surface -> bottom
    ! Arrays are bottom(1) ... surface(N), so we start at N.
    pres(N) = rho(N)*dz(N)/2._rk
    do i = N-1, 1, -1
       pres(i) = pres(i+1) + 0.5_rk * (rho(i+1)*dz(i+1) + rho(i)*dz(i))
    end do

    ! Convert to pressure [dbar] and add atmospheric
    pres(1:N) = p0 + pres(1:N)*gravity/1.0e4_rk

  end subroutine compute_pressure

end module pressure
