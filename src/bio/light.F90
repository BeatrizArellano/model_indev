module light
  use precision_types, only: rk
  implicit none
  private
  public :: compute_PAR_profile, compute_SW_profile

contains

  !--------------------------------------------------------------------
  ! Compute layer-mean downwelling PAR [W m-2] from surface SW (rsds).
  !
  ! Arrays are assumed to describe a *water column only*:
  !   - Local indices 1..N (1 = bottom of this water column, N = surface).
  !   - You can pass slices, e.g. dz(i_seabed:Nfull), par_down(i_seabed:Nfull).
  !
  ! Formula (S2P3-style):
  !   PAR_surface = f_par * rsds
  !   k_par(i)    = k0 + k_chl * chl(i)
  !   par_down(i) = layer-mean PAR in cell i [W m-2]
  !
  ! For now, biology effect is off by default:
  !   k_chl = 0 unless chl_abscross is supplied.
  !--------------------------------------------------------------------
  subroutine compute_PAR_profile(dz, rsds, par_down,  &
                                 par_fraction, par_atten, chl, chl_abscross)

    real(rk), intent(in)           :: dz(:)          ! layer thickness, bottom→surface
    real(rk), intent(in)           :: rsds           ! surface SW downwelling [W m-2]
    real(rk), intent(out)          :: par_down(:)    ! layer-mean PAR [W m-2]
    real(rk), intent(in), optional :: par_fraction   ! fraction of SW that is PAR (~0.45)
    real(rk), intent(in), optional :: par_atten      ! background PAR attenuation [m-1]
    real(rk), intent(in), optional :: chl(:)         ! chlorophyll per layer (optional)
    real(rk), intent(in), optional :: chl_abscross   ! PAR abs. per chl unit [m-1 / (chl unit)]

    integer  :: i, N
    real(rk) :: f_par, k0, k_chl
    real(rk) :: rad_surf, d, dz_i, chl_i
    real(rk), parameter :: tiny = 1.0e-10_rk

    N = size(dz)
    if (size(par_down) /= N) stop 'compute_PAR_profile: size mismatch (par_down)'
    if (present(chl)) then
       if (size(chl) /= N) stop 'compute_PAR_profile: size mismatch (chl)'
    end if

    ! ---- Parameters from S2P3  ----
    f_par = 0.45_rk
    if (present(par_fraction)) f_par = par_fraction

    k0 = 0.10_rk                 ! m-1
    if (present(par_atten)) k0 = par_atten

    k_chl = 0.0_rk               ! biology OFF by default
    if (present(chl_abscross)) k_chl = chl_abscross

    par_down = 0.0_rk

    ! ---- PAR just below the surface (top of layer N) ----
    rad_surf = max(0.0_rk, f_par * rsds)

    ! ---- Integrate from surface (N) to bottom (1) ----
    do i = N, 1, -1
       dz_i = dz(i)

       if (present(chl)) then
          chl_i = chl(i)
       else
          chl_i = 0.0_rk
       end if

       d = max(0.0_rk, k0 + k_chl * chl_i)   ! local k_PAR [m-1]

       if (d * dz_i > tiny) then
          ! Layer-mean PAR = (rad_surf/(d*dz)) * (1 - exp(-d*dz))
          par_down(i) = rad_surf * (1.0_rk - exp(-d*dz_i)) / (d*dz_i)
       else
          ! Very weak attenuation: nearly uniform within the layer
          par_down(i) = rad_surf
       end if

       ! PAR at bottom of this layer (top of the next layer down)
       rad_surf = rad_surf * exp(-d * dz_i)
    end do

  end subroutine compute_PAR_profile


  !--------------------------------------------------------------------
  ! Compute layer-mean downwelling shortwave SW [W m-2].
  !
  ! Arrays are assumed to describe a *water column only*:
  !   - 1 = bottom of the water column, N = surface.
  !
  ! Inputs:
  !   dz(:)          : layer thickness [m], bottom→surface
  !   rsds           : surface downwelling shortwave just below surface [W m-2]
  !
  ! Optional inputs:
  !   chl(:)         : chlorophyll per layer (only used if present)
  !   sw_atten       : background SW attenuation k0 [m-1] (default 0.10)
  !   chl_abscross_sw: SW attenuation per chl unit [m-1 / chl-unit] (default 0.0)
  !
  ! Output:
  !   sw_down(:)     : layer-mean SW irradiance [W m-2], bottom→surface
  !--------------------------------------------------------------------
  subroutine compute_SW_profile(dz, rsds, sw_down, &
                                chl, sw_atten, chl_abscross_sw)
    real(rk), intent(in)           :: dz(:)
    real(rk), intent(in)           :: rsds
    real(rk), intent(out)          :: sw_down(:)
    real(rk), intent(in), optional :: chl(:)
    real(rk), intent(in), optional :: sw_atten, chl_abscross_sw

    integer  :: i, N
    real(rk) :: k0, k_chl
    real(rk) :: rad_surf, d, dz_i, chl_i
    real(rk), parameter :: tiny = 1.0e-10_rk

    N = size(dz)
    if (size(sw_down) /= N) stop 'compute_SW_profile: size mismatch (sw_down)'
    if (present(chl)) then
       if (size(chl) /= N) stop 'compute_SW_profile: size mismatch (chl)'
    end if

    ! ---- S2P3 defaults ----
    k0 = 0.10_rk                  ! background SW attenuation [m-1]
    if (present(sw_atten)) k0 = sw_atten

    k_chl = 0.0_rk                ! biology OFF by default
    if (present(chl_abscross_sw)) k_chl = chl_abscross_sw

    sw_down = 0.0_rk

    ! ---- SW just below the surface (top of layer N) ----
    rad_surf = max(0.0_rk, rsds)

    ! ---- Integrate from surface (N) to bottom (1) ----
    do i = N, 1, -1
       dz_i = dz(i)

       if (present(chl)) then
          chl_i = chl(i)
       else
          chl_i = 0.0_rk
       end if

       d = max(0.0_rk, k0 + k_chl * chl_i)   ! local k_SW [m-1]

       if (d * dz_i > tiny) then
          ! Layer-mean SW
          sw_down(i) = rad_surf * (1.0_rk - exp(-d*dz_i)) / (d*dz_i)
       else
          sw_down(i) = rad_surf
       end if

       ! SW at bottom of this layer (top of the next layer down)
       rad_surf = rad_surf * exp(-d * dz_i)
    end do

  end subroutine compute_SW_profile

end module light
