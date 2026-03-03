!=======================================================================================
!   Shelf Sea Physics:
!   1-D MODEL OF THE EQUATION OF MOTION USING THE Canuto k-e TURBULENCE CLOSURE SCHEME
!
!   This code is fully based on the S2P3 model lineage:
!
!   • Original S2P3 (v7.0):
!       Jonathan Sharples – Univ. of Liverpool & NOC
!       Documented in: Simpson & Sharples (2012), CUP.
!
!   • Regional framework S2P3-R (v1.0):
!       Marsh, Hickman & Sharples (2015), GMD 8, 3163–3178.
!
!   • Large-scale/efficient S2P3-R v2.0:
!       Halloran et al. (2021), GMD 14, 6177–6195.
!
!=======================================================================================
module heat_fluxes
  use forcing_manager, only: ForcingSnapshot
  use physics_params,  only: rho0, rho_air, cp_sw, cp_air, sigma_SB
  use precision_types, only: rk  
  use radiation,       only: compute_radiation_profile

  implicit none
  private
  
  public :: compute_heat_tendency, apply_temperature_tendency

contains

  !---------------------------------------------------------------------------
  ! Compute explicit temperature tendency from surface heat + shortwave penetration
  !
  ! Outputs:
  !   dTdt_heat(1:N) : [K/s] explicit heating tendency for each layer (1=bottom, N=surface)
  !
  ! Optional diagnostics (all [W/m2], positive INTO ocean):
  !   Q_net_surf : net heat flux into the surface cell (nonSW + sw_absorbed_in_top)
  !   max_abs_dTdt : max over layers of abs(dTdt_heat)
  !---------------------------------------------------------------------------
  subroutine compute_heat_tendency(temp, dz, N,                                 &
                                   rsds, rlds_down, wind_speed, airT, rh, airP, &
                                   lw_skin_penetration,                         &
                                   nonvisible_fraction, depth_nonvisible, depth_visible, &
                                   deposit_bottom_residual, apply_bioshading_to_heat, &
                                   dTdt_heat, Q_net_surf, max_abs_dTdt,                    &
                                   swr_c, atten_bio)

      real(rk), intent(in)  :: temp(:)
      real(rk), intent(in)  :: dz(:)
      integer,  intent(in)  :: N
      real(rk), intent(in)  :: rsds, rlds_down, wind_speed, airT, rh, airP
      real(rk), intent(in)  :: lw_skin_penetration, nonvisible_fraction, depth_nonvisible, depth_visible

      logical, intent(in)   :: deposit_bottom_residual, apply_bioshading_to_heat

      real(rk), intent(out)   :: dTdt_heat(1:N)

      real(rk), intent(out), optional :: Q_net_surf
      real(rk), intent(out), optional :: max_abs_dTdt
      real(rk), intent(out)           :: swr_c(1:N)     

      real(rk), intent(in), optional :: atten_bio(1:N)     ! m-1
      
      ! Parameters 
      real(rk), parameter :: eps_sfc = 0.98_rk
      real(rk), parameter :: Cd_h    = 1.45e-3_rk
      real(rk), parameter :: Ce_e    = 1.50e-3_rk

      real(rk) :: atten_local(1:N)
      real(rk) :: swr_abs(1:N)

      real(rk) :: surf_temp
      real(rk) :: sat_vap, vap, spec_hum1, spec_hum2
      real(rk) :: q_up_lw, q_down_lw_eff
      real(rk) :: q_lw_out, q_sensible_out, q_latent_out
      real(rk) :: q_nonSW_in
      real(rk) :: inv_rho_cp
      real(rk) :: maxrate

      atten_local = 0.0_rk
      if (apply_bioshading_to_heat .and. present(atten_bio)) atten_local = atten_bio


      ! Initialise tendency arrays
      dTdt_heat(:) = 0.0_rk
      maxrate = 0.0_rk

      ! Computing shortwave radiation profile
      call compute_radiation_profile(N=N, dz=dz, rsds=rsds, &
                                     nonvisible_fraction=nonvisible_fraction, &
                                     depth_nonvisible=depth_nonvisible, depth_visible=depth_visible, &
                                     apply_bioshading_to_heat=apply_bioshading_to_heat, &
                                     atten_bio=atten_local, &
                                     deposit_bottom_residual=deposit_bottom_residual, &
                                     swr_abs=swr_abs, &
                                     swr_c=swr_c)

      ! --- Surface temperature in Kelvin for Stefan–Boltzmann 
      surf_temp = temp(N) + 273.15_rk

      ! --- Moisture terms ---
      sat_vap = 10.0_rk**((0.7859_rk + 0.03477_rk*temp(N)) / (1.0_rk + 0.00412_rk*temp(N)))  ! mbar
      vap     = max(0.0_rk, 0.01_rk*rh*sat_vap)                                              ! mbar

      spec_hum1 = 0.62_rk*sat_vap/(airP - 0.38_rk*sat_vap)  ! specific humidity at surface 
      spec_hum2 = 0.62_rk*vap/(airP - 0.38_rk*vap)          ! specific humidity of air 

      ! --- Longwave, sensible, latent: positive out of the sea (S2P3)  ---
      q_up_lw       = eps_sfc * sigma_SB * surf_temp**4
      q_down_lw_eff = lw_skin_penetration * rlds_down
      q_lw_out      = q_up_lw - q_down_lw_eff

      ! --- Sensible/latent (positive OUT of sea) ---
      q_sensible_out = Cd_h * rho_air * cp_air * wind_speed * (temp(N) - airT)
      q_latent_out   = Ce_e * rho_air * wind_speed * (spec_hum1 - spec_hum2) * (2.5e6_rk - 2.3e3_rk*temp(N))

      ! Convert to "positive INTO ocean"
      q_nonSW_in = -(q_lw_out + q_sensible_out + q_latent_out)   

      inv_rho_cp = 1.0_rk / (rho0 * cp_sw)

      ! --- Surface cell tendency: (nonSW + sw absorbed in top) / (rho*cp*dz) ---
      ! SWR absorption heats all layers
      dTdt_heat(1:N) = dTdt_heat(1:N) + (swr_abs(1:N) * inv_rho_cp) / dz(1:N)

      ! Non-shortwave fluxes act at the surface cell only
      dTdt_heat(N) = dTdt_heat(N) + (q_nonSW_in * inv_rho_cp) / dz(N)

      if (present(Q_net_surf)) Q_net_surf = q_nonSW_in + swr_abs(N)

      maxrate = maxval(abs(dTdt_heat(1:N)))
      if (present(max_abs_dTdt)) max_abs_dTdt = maxrate

  end subroutine compute_heat_tendency


  subroutine apply_temperature_tendency(temp, N, dt, dTdt_heat)
      real(rk), intent(inout) :: temp(1:N)
      integer,  intent(in)    :: N
      real(rk), intent(in)    :: dt
      real(rk), intent(in)    :: dTdt_heat(1:N)

      temp(1:N) = temp(1:N) + dt * dTdt_heat(1:N)
  end subroutine apply_temperature_tendency


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
