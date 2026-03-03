!===============================================================================
! radiation.F90
!
! Two-band shortwave radiation scheme for physics:
!   - Non-visible band (UV+IR) with shallow e-folding depth
!   - Visible band with deeper e-folding depth + optional extra attenuation
!     from bioshading (derived from FABM and passed in as atten_bio)
!
! Produces:
!   swr_dn(0:N)   : downward SWR flux at interfaces [W m-2]
!   swr_abs(1:N)  : SWR absorbed in each layer      [W m-2]
!   swr_c(1:N)    : SWR at layer centres (diagnostic) [W m-2]
!
! Index convention:
!   Layers are 1..N with 1=bottom, N=surface.
!   Interfaces are 0..N with N=surface interface, 0=bottom interface.
!===============================================================================
module radiation
  use precision_types, only: rk
  implicit none
  private

  public :: compute_radiation_profile, compute_par_profile

contains

    !---------------------------------------------------------------------------
    ! Compute two-band shortwave radiation profiles.
    !
    ! Inputs:
    !   N, dz(1:N)            : layer thickness [m], bottom->surface
    !   rsds                  : downwelling SWR at surface [W m-2]
    !
    !   nonvisible_fraction   : fraction of rsds in non-visible band (0..1)
    !   depth_nonvisible      : e-folding depth of non-visible band [m]
    !   depth_visible         : background e-folding depth of visible band [m]
    !
    ! Optional bioshading for HEAT:
    !   apply_bioshading_to_heat : if true, visible attenuation gets +atten_bio(i)
    !   atten_bio(1:N)          : extra attenuation [m-1] (e.g. from FABM)
    !
    ! Bottom handling:
    !   deposit_bottom_residual : if true, add swr_dn(0) to swr_abs(1) and set swr_dn(0)=0
    !
    ! Outputs:
    !   swr_dn(0:N)  : downward SWR flux at interfaces [W m-2]
    !   swr_abs(1:N) : absorbed SWR in each layer [W m-2]
    !   swr_c(1:N)   : SWR at layer centres [W m-2] (diagnostic)
    !
    ! Notes:
    !   - Energy conservation (water column):
    !       sum_i swr_abs(i) + swr_dn(0) = rsds
    !     If deposit_bottom_residual=.true. then sum_i swr_abs(i) = rsds.
    !
    !   - swr_c is computed as centre value of total SWR (visible + non-visible):
    !       swr_c(i) = F_top(i) * exp(-k_eff(i) * 0.5*dz(i))
    !     where k_eff is a flux-weighted combination of band attenuations.
    !     This avoids needing to output band-separated centre values.
    !---------------------------------------------------------------------------
    subroutine compute_radiation_profile(N, dz, rsds,                                         &
                                        nonvisible_fraction, depth_nonvisible, depth_visible, &
                                        deposit_bottom_residual,                              &
                                        swr_abs, swr_c, apply_bioshading_to_heat, atten_bio)

        integer,  intent(in)  :: N
        real(rk), intent(in)  :: dz(1:N)
        real(rk), intent(in)  :: rsds

        real(rk), intent(in)  :: nonvisible_fraction
        real(rk), intent(in)  :: depth_nonvisible
        real(rk), intent(in)  :: depth_visible

        logical,  intent(in)  :: deposit_bottom_residual

        real(rk), intent(out) :: swr_abs(1:N)
        real(rk), intent(out) :: swr_c(1:N)

        logical,  intent(in)  :: apply_bioshading_to_heat
        real(rk), intent(in), optional :: atten_bio(1:N)

        integer  :: i
        real(rk) :: f_nv, f_vis
        real(rk) :: k_nv, k_vis_bg, k_vis
        real(rk) :: Fin_nv, Fout_nv, Fin_vis, Fout_vis
        real(rk) :: Ftop_total
        real(rk) :: kbio, k_eff
        real(rk) :: residual
        real(rk), parameter :: tiny = 1.0e-12_rk


        f_nv  = nonvisible_fraction
        f_vis = 1.0_rk - f_nv

        k_nv     = 1.0_rk / max(depth_nonvisible, 1.0e-6_rk)
        k_vis_bg = 1.0_rk / max(depth_visible,    1.0e-6_rk)

        ! initialise outputs
        swr_abs(:) = 0.0_rk
        swr_c(:)   = 0.0_rk

        ! band-separated fluxes at top
        Fin_nv  = f_nv  * max(0.0_rk, rsds)
        Fin_vis = f_vis * max(0.0_rk, rsds)

        ! integrate from surface (N) to bottom (1)
        do i = N, 1, -1

            kbio = 0.0_rk
            if (apply_bioshading_to_heat .and. present(atten_bio)) then
                kbio = max(0.0_rk, atten_bio(i))
            end if
            k_vis = max(0.0_rk, k_vis_bg + kbio)

            ! propagate each band across layer i
            Fout_nv  = Fin_nv  * exp(-k_nv  * dz(i))
            Fout_vis = Fin_vis * exp(-k_vis * dz(i))

            ! absorbed energy in layer i (W/m2)
            swr_abs(i) = (Fin_nv - Fout_nv) + (Fin_vis - Fout_vis)


            ! Computing the radiation profile at layers' centres
            Ftop_total = Fin_nv + Fin_vis
            if (Ftop_total > tiny) then
                k_eff = max(0.0_rk, (Fin_nv*k_nv + Fin_vis*k_vis) / Ftop_total)
                swr_c(i) = Ftop_total * exp(-k_eff * 0.5_rk * dz(i))
            else
                swr_c(i) = 0.0_rk
            end if

            ! advance downward
            Fin_nv  = Fout_nv
            Fin_vis = Fout_vis

        end do

        ! Optionally deposit any residual reaching the seabed into the bottom layer
        residual = Fin_nv + Fin_vis
        if (deposit_bottom_residual) swr_abs(1) = swr_abs(1) + residual

    end subroutine compute_radiation_profile

    subroutine compute_par_profile(N, dz, rsds,                          &
                                   par_fraction, depth_visible,          &
                                   apply_bioshading_to_par, atten_bio,   &
                                   par, par_sfc)

        integer,  intent(in)  :: N
        real(rk), intent(in)  :: dz(1:N)
        real(rk), intent(in)  :: rsds

        real(rk), intent(in)  :: par_fraction 
        real(rk), intent(in)  :: depth_visible           ! Background PAR e-folding depth [m]

        logical,  intent(in)  :: apply_bioshading_to_par
        real(rk), intent(in), optional :: atten_bio(1:N) ! Attenuation from organisms/particles [m-1]

        real(rk), intent(out) :: par(1:N)          ! PAR at layer centres [W m-2]
        real(rk), intent(out), optional :: par_sfc ! PAR at the surface

        integer :: i
        real(rk) :: k_par_bg, kbio, k_par
        real(rk) :: Fin_par, Fout_par
        real(rk), parameter :: tiny = 1.0e-12_rk

        if (par_fraction < 0.0_rk .or. par_fraction > 1.0_rk) then
            stop 'compute_par_profile: par_fraction must be in [0,1]'
        end if

        k_par_bg = 1.0_rk / max(depth_visible, 1.0e-6_rk)

        par(:) = 0.0_rk

        ! Downward PAR flux just below surface (top of layer N)
        Fin_par = par_fraction * max(0.0_rk, rsds)
        if (present(par_sfc)) par_sfc = Fin_par

        do i = N, 1, -1

            kbio = 0.0_rk
            if (apply_bioshading_to_par .and. present(atten_bio)) then
                kbio = max(0.0_rk, atten_bio(i))
            end if
            k_par = max(0.0_rk, k_par_bg + kbio)

            ! Centre value
            if (Fin_par > tiny) then
                par(i) = Fin_par * exp(-k_par * 0.5_rk * dz(i))
            else
                par(i) = 0.0_rk
            end if

            ! propagate to bottom of layer for next step down
            Fout_par = Fin_par * exp(-k_par * dz(i))
            Fin_par  = Fout_par

        end do

    end subroutine compute_par_profile

end module radiation