module heat_fluxes
  use precision_types, only: rk
  use physics_params,  only: rho0, rho_air, cp_sw, cp_air, sigma_SB
  use forcing_manager, only: ForcingSnapshot
  implicit none
  private
  public :: SURFACE_HEAT

contains

  !---------------------------------------------------------------------------
  ! SURFACE_HEAT
  !
  ! Sign convention: positive Q means heat INTO the ocean.
  !
  ! Inputs (per time step):
  !   temp_old(N)     : old potential temperature [°C]
  !   dz(N)           : layer thickness [m]
  !   dt              : time step [s]
  !   rsds            : shortwave downward flux at surface (0–) [W m-2]
  !   rlds_down       : longwave downward atmospheric flux (0–) [W m-2]
  !   wind_speed      : 10 m wind speed [m s-1]
  !   airT            : 2 m air temperature [°C]
  !   rh              : relative humidity [%] (0–100)
  !   airP            : sea-level pressure [hPa]
  !   lambda          : base light attenuation coefficient [m-1]
  !   heat_shade      : attenuation per (chl) unit [m-1 per chl-unit]
  !   chl(N)          : proxy for x_new(i) in S2P3 (e.g., chlorophyll) [same units your model uses]
  ! lw_skin_penetration: factor applied to downwelling LW to mimic skin absorption (S2P3 used 0.95)
  !
  !---------------------------------------------------------------------------
  subroutine SURFACE_HEAT(temp, dz, N, dt,                         &
                          rsds, rlds_down, wind_speed, airT, rh, airP, &
                          lw_skin_penetration,lambda,  &
                          heat_shade, chl)

    real(rk), intent(inout) :: temp(:)
    real(rk), intent(in)    :: dz(:)
    integer,  intent(in)    :: N
    real(rk), intent(in)    :: dt
    real(rk), intent(in)    :: rsds, rlds_down, wind_speed, airT, rh, airP
    real(rk), intent(in)    :: lambda, lw_skin_penetration
    real(rk), intent(in), optional :: heat_shade
    real(rk), intent(in), optional :: chl(:)

    !real(rk), intent(out) :: temp_new(:)
    

    real(rk), parameter :: eps_sfc   = 0.98_rk   ! Ocean surface emissivity for longwave radiation (fraction)
    real(rk), parameter :: Cd_h      = 1.45e-3_rk
    real(rk), parameter :: Ce_e      = 1.50e-3_rk
    real(rk), parameter :: sw_top_fraction = 0.55_rk
    ! -------------------------------------------------------------------------

    integer  :: i
    real(rk) :: surf_temp, sat_vap, vap, spec_hum2, spec_hum1
    real(rk) :: q_up_lw, q_down_lw_eff
    real(rk) :: q_sensible, q_latent
    real(rk) :: q_sw_top, q_sw_rest, flux, atten, Iin, Iout
    real(rk) :: q_lw, rad_out, rad_net, inv_rho_cp_dt
    real(rk) :: k_top 
    real(rk) :: chl_i, shade

    !temp_new(:) = temp_old(:)

    if (present(heat_shade)) then
      shade = heat_shade
    else
      shade = 0.012_rk
    end if
    if (present(chl)) then
      chl_i = chl(N)
    else
      chl_i = 0.0_rk
    end if

    surf_temp = (temp(N) + 273.15_rk)  ! Skin temperatyre in °K

    ! --- Moisture terms --- 
    sat_vap = 10.0_rk**((0.7859_rk + 0.03477_rk*temp(N)) / (1.0_rk + 0.00412_rk*temp(N)))   ! mbar
    vap     = max(0.0_rk, 0.01_rk*rh*sat_vap)                                                       ! mbar
    ! specific humidity approximations
    spec_hum1 = 0.62_rk*sat_vap/(airP - 0.38_rk*sat_vap)  ! specific humidity at surface 
    spec_hum2 = 0.62_rk*vap/(airP - 0.38_rk*vap)          ! specific humidity of air 
    ! --- Longwave, sensible, latent: positive out of the sea (S2P3)  ---
    q_up_lw       = eps_sfc * sigma_SB * surf_temp**4
    q_down_lw_eff = lw_skin_penetration * rlds_down
    q_lw          = q_up_lw - q_down_lw_eff

    ! Sensible and latent 
    q_sensible = Cd_h * rho_air * cp_air * wind_speed * (temp(N) - airT)          
    q_latent   = Ce_e * rho_air * wind_speed * (spec_hum1 - spec_hum2) * (2.5e6_rk - 2.3e3_rk*temp(N))

    rad_out = -(q_lw + q_sensible + q_latent)

    ! --- Shortwave split: sw_top_fraction to top layer, remainder decays with Beer–Lambert ---
    !q_sw_top  = sw_top_fraction * rsds               ! Legacy S2P3
    !q_sw_rest = (1.0_rk - sw_top_fraction) * rsds

    ! Trying here a physical-based absorption in the top instead
    k_top = lambda + shade*chl_i
    ! Apply top-cell shortwave directly (as heat tendency)
    q_sw_top = rsds * (1 - exp(-k_top * dz(N)))
    q_sw_rest = rsds - q_sw_top   
  
    inv_rho_cp_dt = dt / (rho0 * cp_sw)    

    ! Update surface temperature
    temp(N) = temp(N)+ (rad_out+q_sw_top)*inv_rho_cp_dt/dz(N)     ! Radiation lost at the surface plus shortwave there

    ! Distribute remaining shortwave through the column from surface downward decaying with Beer-Lambert
    Iin = q_sw_rest

    do i = N-1, 1, -1
      if (present(chl)) chl_i = chl(i)  ! <-- use chl if provided, else 0
      atten   = max(0.0_rk, lambda + shade * chl_i)   ! local attenuation
      Iout = Iin * exp(-atten * dz(i))
      flux = (Iin - Iout)                                ! W m-2 absorbed in layer i
      temp(i) = temp(i) + inv_rho_cp_dt * (flux / dz(i))
      Iin  = Iout
    end do
    ! Any residue below bottom goes to the bottom layer
    if (Iin > 0.0_rk) then
      temp(1) = temp(1) + inv_rho_cp_dt * (Iin / dz(1))
      Iin = 0.0_rk
    end if

  end subroutine SURFACE_HEAT


end module heat_fluxes
