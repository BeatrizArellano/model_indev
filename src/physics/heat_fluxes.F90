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
  use physics_params,  only: rho0, cp_sw, cp_air, sigma_SB 
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
  subroutine compute_heat_tendency(temp, dz, N,                                          &
                                   rsds, rlds_down, wind_speed, airT, rh, airP,          &
                                   nonvisible_fraction, depth_nonvisible, depth_visible, &
                                   deposit_bottom_residual, apply_bioshading_to_heat,    &
                                   dTdt_heat, Q_net_surf, Q_nonSW_in, max_abs_dTdt,      &
                                   swr_c, atten_bio)

      real(rk), intent(in)  :: temp(:)
      real(rk), intent(in)  :: dz(:)
      integer,  intent(in)  :: N
      real(rk), intent(in)  :: rsds, rlds_down, wind_speed, airT, rh, airP         ! airP : sea-level pressure [hPa]
      real(rk), intent(in)  :: nonvisible_fraction, depth_nonvisible, depth_visible

      logical, intent(in)   :: deposit_bottom_residual, apply_bioshading_to_heat

      real(rk), intent(out)   :: dTdt_heat(1:N)

      real(rk), intent(out), optional :: Q_net_surf
      real(rk), intent(out), optional :: Q_nonSW_in
      real(rk), intent(out), optional :: max_abs_dTdt
      real(rk), intent(out)           :: swr_c(1:N)     

      real(rk), intent(in), optional :: atten_bio(1:N)     ! m-1
      
      ! Parameters 
      real(rk), parameter :: emiss_sea = 0.985_rk

      real(rk) :: atten_local(1:N)
      real(rk) :: swr_abs(1:N)

      real(rk) :: surf_temp
      real(rk) :: sat_vap_air, sat_vap_sea, vap_air, spec_hum1, spec_hum2
      real(rk) :: q_up_lw
      real(rk) :: q_lw_out, q_sensible_out, q_latent_out
      real(rk) :: q_nonSW_in_local
      real(rk) :: Ch_h, Ce_q
      real(rk) :: inv_rho_cp, rho_air
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
      ! Arden Buck empirical fit
      sat_vap_air = 6.1121_rk * exp((18.678_rk - airT/234.5_rk) * (airT/(257.14_rk + airT)))          ! mbar
      sat_vap_sea = 6.1121_rk * exp((18.678_rk - temp(N)/234.5_rk) * (temp(N)/(257.14_rk + temp(N))))
      vap_air     = max(0.0_rk, 0.01_rk * rh * sat_vap_air)                                           ! mbar

      spec_hum1 = 0.622_rk*sat_vap_sea/(airP - 0.378_rk*sat_vap_sea)  ! specific humidity at surface 
      spec_hum2 = 0.622_rk*vap_air    /(airP - 0.378_rk*vap_air)      ! specific humidity of air 

      ! --- Longwave, sensible, latent: positive out of the sea (S2P3)  ---
      q_up_lw  = emiss_sea * sigma_SB * surf_temp**4
      q_lw_out = q_up_lw - emiss_sea * rlds_down

      ! Dynamic bulk transfer coefficients for moisture and heat
      Ce_q = (1.05e-3_rk + 0.05e-3_rk * wind_speed)    ! bulk transfer coefficient for latent heat
      Ch_h = (0.75e-3_rk + 0.04e-3_rk * wind_speed)    ! bulk transfer coefficient for sensible heat
      
      Ce_q = min(max(Ce_q, 1.05e-3_rk), 1.60e-3_rk)
      Ch_h = min(max(Ch_h, 0.75e-3_rk), 1.25e-3_rk)

      ! Scaling air density with temperature
      rho_air = airP * 100.0_rk / (287.05_rk * (airT + 273.15_rk))

      ! --- Sensible/latent (positive OUT of sea) ---
      q_sensible_out = Ch_h * rho_air * cp_air * wind_speed * (temp(N) - airT)
      q_latent_out   = Ce_q * rho_air * wind_speed * (spec_hum1 - spec_hum2) * (2.5e6_rk - 2.3e3_rk*temp(N))

      ! Convert to "positive INTO ocean"
      q_nonSW_in_local = -(q_lw_out + q_sensible_out + q_latent_out)   

      inv_rho_cp = 1.0_rk / (rho0 * cp_sw)

      ! --- Surface cell tendency: (nonSW + sw absorbed in top) / (rho*cp*dz) ---
      ! SWR absorption heats all layers
      dTdt_heat(1:N) = dTdt_heat(1:N) + (swr_abs(1:N) * inv_rho_cp) / dz(1:N)

      if (present(Q_nonSW_in)) Q_nonSW_in = q_nonSW_in_local
      if (present(Q_net_surf)) Q_net_surf = q_nonSW_in_local + swr_abs(N)

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

end module heat_fluxes
