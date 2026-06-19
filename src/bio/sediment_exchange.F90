module sediment_exchange
    use precision_types,   only: rk
    implicit none
    private

    public :: compute_solute_flux_swi
    public :: apply_particulate_deposition
    public :: apply_bioirrigation

    ! Kinematic viscosity of water (ν) [m^2 s^-1] (used as a Constant)
    ! real(rk), parameter :: kin_visc = 1.3e-6_rk
    
    real(rk), parameter :: kappa = 0.4_rk   ! von Kármán constant [-]

contains

    !---------------------------------------------------------------------------------
    ! Compute hydrodynamically-controlled solute flux across the sediment-water interface
    ! following Umlauf et al. (2023) approach.
    !
    ! Sign convention:
    !   swi_flux > 0  => flux from sediment to water (upward, source to the pelagic bottom layer)
    !   swi_flux < 0  => flux from water into sediment (downward, sink from the pelagic bottom layer)
    !
    ! Inputs:
    !   c_w     : solute concentration in bottom water layer (layer-centre)        [Var]
    !   c_s     : solute concentration in top sediment layer (layer-centre)        [Var]
    !
    !   u_taub  : bottom friction velocity u*                                      [m s^-1]
    !            (sets turbulent shear scale that controls exchange intensity)
    !
    !   z0b     : bottom roughness length z0                                      [m]
    !            (hydrodynamic roughness controlling log-law and transfer)
    !
    !   h0b     : Physical roughness height   (ks)                               [m]
    !
    !   Nz      : turbulent diffusivity of momentum at the top of the 
    !             bottom layer of water column      [m^2 s^-1]
    !
    !   Kz      : turbulent diffusivity for scalars at the top of the 
    !             bottom layer of water column   [m^2 s^-1]
    !
    !   D_sol   : Molecular diffusivity of the solute  [m^2 s^-1]
    !
    !   D_eff   : effective porewater diffusivity in sediment (tortuosity-corrected)  [m^2 s^-1]
    !
    !   dz_sed  : thickness of top sediment layer                                 [m]   
    !
    !   dz_w     : thickness of bottom water layer                                [m]
    !
    !   kin_visc : kinematic viscosity of bottom seawater [m2 s-1]
    !
    ! Outputs:
    !   swi_flux : Flux at the sediment-water interface (positive upward into water) [Var m^-2 s^-1]
    !              This flux is used as a Neumann boundary condition for diffusion at the
    !              bottom of the water column and at the top of the sediment. 
    !
    ! Reference:   Umlauf, L., Klingbeil, K., Radtke, H., Schwefel, R., Bruggeman, J., & Holtermann, P.
    !              (2023). Hydrodynamic Control of Sediment‐Water Fluxes: Consistent Parameterization 
    !              and Impact in Coupled Benthic‐Pelagic Models. Journal of Geophysical Research: Oceans, 
    !              128(6), e2023JC019651.       
    !
    !   Umlauf et al. (2023) replace the traditional approach of prescribing a
    !   diffusive boundary layer (DBL) thickness with a hydrodynamically controlled
    !   exchange formulation. A key motivation for this is to remove artificial 
    !   grid dependence in SWI fluxes. Instead of assuming a fixed DBL, 
    !   the sediment–water flux emerges from the balance between turbulent transport 
    !   at the bottm boundary layer and molecular diffusion in the sediment.
    !   Bottom stress (via the friction velocity) controls the water-side
    !   transport dominated by turbulece, while porewater molecular diffusivity controls 
    !   the sediment-side transport.
    !   The transfer coefficient beta is evaluated using the large-Schmidt-number approximations
    !   of Umlauf et al. (2023). All biogeochemical tracers have large Schmidt numbers. 
    !   The smooth, transitional, and rough regimes are distinguished using the roughness 
    !   Reynolds number ks+ = ks*u*/nu, where ks is the equivalent sand-grain roughness height.
    !---------------------------------------------------------------------------------
    subroutine compute_solute_flux_swi(c_w, c_s, phi_swi, u_taub, z0b, h0b,    &
                                       Nz, Kz, D_sol, D_eff, dz_w, dz_sed,      &
                                       kin_visc, swi_flux)

        ! Inputs
        real(rk), intent(in) :: c_w, c_s
        real(rk), intent(in) :: phi_swi             ! Porosity at the sediment-watre interface
        real(rk), intent(in) :: u_taub, z0b, h0b, dz_w
        real(rk), intent(in) :: Nz, Kz
        real(rk), intent(in) :: D_sol, D_eff
        real(rk), intent(in) :: dz_sed
        real(rk), intent(in) :: kin_visc            ! Kinematic viscosity of water (ν) [m^2 s^-1]

        ! Outputs
        real(rk), intent(out) :: swi_flux       

        ! Locals
        real(rk) :: L_wat, L_sed, logterm
        real(rk) :: Pr, Sc, beta, r_c, G_wat, G_sed, c_swi
        real(rk) :: ks_plus, z0r_plus, beta_smooth, beta_rough, alpha

        real(rk), parameter :: B_smooth = 5.5_rk
        real(rk), parameter :: B_rough  = 8.5_rk
        real(rk), parameter :: eps = 1.0e-20_rk

        !--------------------------
        ! Turbulent Prandtl number Pr = K_m / K_h  (momentum vs scalar mixing)
        ! at the top of the bottom layer of the water column
        ! From Eq. 13 in Umlauf et al. (2023)
        !--------------------------
        Pr = Nz / max(Kz, eps)

        !--------------------------
        ! Molecular Schmidt/Prandtl number for the solute: Sc = nu / D
        !  (ratio of momentum diffusivity to molecular diffusivity)
        ! Appendix B in Umlauf et al. (2023)
        !--------------------------
        Sc = kin_visc / max(D_sol, eps)

        ! Roughness Reynolds number based on physical roughness height
        ! Yaglom & Kader use ks+ = ks*u*/nu for the transitional interpolation.
        ks_plus = h0b * u_taub / max(kin_visc, eps)
        z0r_plus = (h0b / 30.0_rk) * u_taub / max(kin_visc, eps)

        !--------------------------
        ! Beta function
        ! Approximation for Smooth-bottom only
        ! Appendix B in Umlauf et al. (2023), eq B2
        !--------------------------
        beta_smooth = 14.8_rk * Sc**(2.0_rk/3.0_rk)        

        !-------------------------------------------
        ! Beta function
        ! Approximation for rough-bottom only
        ! Appendix B in Umlauf et al. (2023), eq B3
        !-------------------------------------------
        beta_rough = 0.55_rk * exp(kappa * B_rough / 2.0_rk) * &
                     sqrt(z0r_plus) * Sc**(2.0_rk/3.0_rk) !- Pr * B_rough + 9.5_rk

        ! Smooth regime
        if (ks_plus <= 3.0_rk) then           ! Consistent with z0+ < 0.1
            beta = beta_smooth
        ! Rough bottom regime
        else if (ks_plus >= 100.0_rk) then    ! Fully rough following Yaglom & Kader; z0r+ ≈ 3.3
            beta = beta_rough
        ! Transitional regime (3 > ks+ < 100)
        else
            ! Transitional interpolation adapted from Yaglom & Kader (1974) eq. 23
            alpha = (ks_plus - 3.0_rk) / (100.0_rk - 3.0_rk)
            beta = (1.0_rk - alpha) * beta_smooth + alpha * beta_rough
        end if

        !--------------------------
        ! Transfer function r_c(L_wat):
        ! Eq. 17 in Umlauf et al. (2023)
        ! rc(L_wat) = 1/[(Sc_t/kappa)*ln((L_wat+z0)/z0) + beta]
        !--------------------------
        L_wat = 0.5_rk * dz_w              ! Height above seabed: centre of bottom water layer
        logterm = log(L_wat / max(z0b, eps) + 1.0_rk )
        r_c     = 1.0_rk / max((Pr / kappa) * logterm + beta, eps)   ! dimensionless

        !--------------------------
        ! Water-side diffusive conductance (transfer velocity) G_wat  [m/s]
        ! From Eq. 16 in Umlauf et al. (2023)
        ! How efficiently turbulence can transport a solute between the SWI and the bottom water layer. 
        !--------------------------
        G_wat = u_taub * r_c

        !--------------------------
        ! Sediment-side diffusive conductance of the top sediment layer G_sed [m/s]
        ! From Eq. 22 in Umlauf et al. (2023)
        !   G_sed ~ D_eff / (0.5*dz_sed)
        !--------------------------
        L_sed = 0.5_rk * dz_sed     ! Distance between the centre of the top sediment layer and the sediment-water interface
        G_sed = max(phi_swi, eps) * D_eff / L_sed

        !--------------------------
        ! Concentration at the sediment-water interface
        ! Eq. 23 in Umlauf et al. (2023)
        !  c_swi = (G_sed*c_s + G_wat*c_w) / (G_sed + G_wat)
        !--------------------------
        c_swi = (G_sed*c_s + G_wat*c_w) / max(G_sed + G_wat, eps)

        !--------------------------
        ! Flux at the sediment-water interface (positive upward into water)
        ! Eq. 16 in Umlauf et al. (2023)
        !--------------------------
        swi_flux = G_wat * (c_swi - c_w)

    end subroutine compute_solute_flux_swi


    !! Apply particulate deposition across the sediment–water interface (SWI).
    !!
    !! Transfers particulate material from the bottom water cell into the top
    !! sediment layer based on the downward settling velocity at the SWI.
    !! Mass is conserved by removing material from the water column (per water
    !! volume) and adding it to the sediment inventory expressed per bulk
    !! sediment volume. Only downward (negative) velocities contribute to
    !! deposition; no resuspension or upward exchange is represented.
    subroutine apply_particulate_deposition(Cw_bot, dz_w_bot, vel_swi, dt, dz_sed_top, Cbulk_sed_top)
        real(rk), intent(inout) :: Cw_bot          ! bottom-water concentration (per water volume)
        real(rk), intent(in)    :: dz_w_bot        ! thickness of bottom-water cell [m]
        real(rk), intent(in)    :: vel_swi         ! velocity at SWI interface [m/s] (negative downward)
        real(rk), intent(in)    :: dt              ! [s]
        real(rk), intent(in)    :: dz_sed_top      ! thickness of top sediment layer [m]
        real(rk), intent(inout) :: Cbulk_sed_top   ! bulk-sediment conc in top layer (per bulk volume)

        real(rk) :: dep_flux_req, dep_mass_req, dep_mass
        real(rk) :: available

        ! Only downward movement deposits into sediment.
        if (vel_swi < 0._rk) then
            ! Available particulate inventory in bottom water cell [mass m-2].
            available = max(Cw_bot, 0._rk) * dz_w_bot
            if (available <= 0._rk) return

            ! Requested depositional flux [mass m-2 s-1].
            dep_flux_req = (-vel_swi) * max(Cw_bot, 0.0_rk)

            ! Requested depositional inventory over this time step [mass m-2].
            dep_mass_req = dt * dep_flux_req

            ! Positivity-preserving deposited inventory.
            dep_mass = min(dep_mass_req, available)

            ! Update bottom-water concentration.
            Cw_bot = Cw_bot - dep_mass / dz_w_bot

            ! Avoid tiny negative roundoff.
            if (Cw_bot < 0.0_rk .and. Cw_bot > - 100._rk * epsilon(1._rk)) Cw_bot = 0.0_rk

            ! update sediment surface concentration (Forward Euler)
            Cbulk_sed_top = Cbulk_sed_top + dep_mass / dz_sed_top
        end if
    end subroutine apply_particulate_deposition

    !-----------------------------------------------------------------------
    ! apply_bioirrigation: Apply bioirrigation as a non-local exchange 
    ! between sediment porewater and the overlying bottom-water layer.
    !
    ! Concept:
    !   Canonical early-diagenesis models often represent bioirrigation as
    !   a first-order volumetric exchange of porewater with bottom water:
    !
    !     dC_s(k)/dt = alpha(k) * (C_bw - C_s(k))
    !
    !   where alpha(k) [s^-1] is the fraction of porewater volume in layer k
    !   exchanged per unit time (typically decaying with depth). This is a
    !   standard “non-local exchange” closure used to represent burrow
    !   flushing/ventilation without resolving burrow geometry.
    !
    ! Discrete, mass-conserving interpretation (per unit horizontal area):
    !   Porewater volume per unit area in layer k is
    !     Hpw(k) = phi(k) * dz_sed(k)  [m]  (stored as porewat_thickness)
    !
    !   Over a time step dt, exchanged porewater volume per area is:
    !     dV_k = alpha(k) * Hpw(k) * dt
    !
    !   Net mass gained by sediment layer k (and lost by bottom water) is:
    !     dM_k = (C_bw - C_s(k)) * dV_k
    !          = alpha(k) * (C_bw - C_s(k)) * Hpw(k) * dt   [tracer units m s-1]
    !
    !   Updating sediment and bottom water with equal-and-opposite dM_k
    !   ensures exact tracer mass conservation (up to roundoff).
    !
    ! Notes:
    !   - concentration(1:nsed) must be porewater concentrations for solutes.
    !   - concentration(k_wat_btm) is the bottom-water concentration in the
    !     first water layer above the SWI.
    !   - This explicit update is conservative but not positivity-preserving;
    !     checks below require sufficiently small dt (or global substepping).
    ! 
    ! Based on section 3.10.5 in Boudreau (1997) and implementation of
    ! equation 3.140 therein to model bioirrigation as a non-local exchange. 
    ! Specific discussion of this is in Boudreau (1984)
    !-----------------------------------------------------------------------
    subroutine apply_bioirrigation(dt, nsed, ntotal, alpha, porewat_thickness, dz_wat_btm, k_wat_btm,   &
                                   concentration, bioirr_flux)
        real(rk), intent(in)    :: dt
        integer,  intent(in)    :: nsed, ntotal
        real(rk), intent(in)    :: alpha(nsed)               ! Bioirrigation alpha [s^-1]       
        real(rk), intent(in)    :: porewat_thickness(nsed)   ! [m] porewater capacity per m2
        real(rk), intent(in)    :: dz_wat_btm                ! [m] Thickness of the bottom layer in water column        
        integer,  intent(in)    :: k_wat_btm                 ! Indices for water bottom and sediment 
        real(rk), intent(inout) :: concentration(ntotal)     ! Concentration per porewater volume in sediments and water volume in the water column        
        real(rk), intent(out)   :: bioirr_flux               ! Sediment-water flux due to bioirrigation [concentration units * m s-1]

        integer  :: k
        real(rk) :: Hpw, alpha_k
        real(rk) :: Cbw, Cs, dM, dMtot        

        Cbw = concentration(k_wat_btm)
        dMtot = 0._rk
        bioirr_flux = 0._rk

        do k = 1, nsed
            alpha_k = alpha(k)
            if (alpha_k <= 0._rk) cycle

            ! The amount of mass that can be exchanged depends on how much porewater exists in that layer.
            ! Hpw is the actual porewater volume available to be exchanged
            Hpw = porewat_thickness(k)                 ! [m] porewater capacity per m2
            if (Hpw <= 0._rk) cycle

            Cs = concentration(k)                      ! porewater concentration 

            ! Compute how much mass per m2 is exchanged over dt for each sediment layer
            ! The exchanged volume per unit area is ΔVk = alpha_k* Hpw * dt 
            !   ΔM_k = (Cbw - Ck) * ΔVk
            ! Exchange mass per unit area [concentration units * m] 
            dM = alpha_k * (Cbw - Cs) * Hpw * dt

            ! Update sediment porewater concentration in the kth sediment layer
            ! by distributing the mass over its porewater volume per m2
            concentration(k) = Cs + dM / Hpw

            ! Accumulate mass removed (or added) from bottom water (opposite sign)
            dMtot = dMtot - dM

            ! Sediment-water flux due to bioirrigation [concentration units * m s-1]
            ! Positive into the water, negative into the sediments
            bioirr_flux = bioirr_flux - dM / dt
        end do

        ! Update bottom water concentration
        concentration(k_wat_btm) = Cbw + dMtot / dz_wat_btm

    end subroutine apply_bioirrigation

end module sediment_exchange