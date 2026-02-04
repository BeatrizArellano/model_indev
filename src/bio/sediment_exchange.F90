module sediment_exchange
    use precision_types,   only: rk
    implicit none
    private

    public :: compute_solute_flux_swi
    public :: apply_particulate_deposition
    public :: apply_bioirrigation

    ! Kinematic viscosity of water (ν) [m^2 s^-1] (used as a Constant)
    real(rk), parameter :: kin_visc = 1.3e-6_rk
    
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
    !---------------------------------------------------------------------------------
    subroutine compute_solute_flux_swi(c_w, c_s, phi_swi, u_taub, z0b,         &
                                       Nz, Kz, D_sol, D_eff, dz_w, dz_sed,   &
                                       swi_flux)

        ! Inputs
        real(rk), intent(in) :: c_w, c_s
        real(rk), intent(in) :: phi_swi             ! Porosity at the sediment-watre interface
        real(rk), intent(in) :: u_taub, z0b, dz_w
        real(rk), intent(in) :: Nz, Kz
        real(rk), intent(in) :: D_sol, D_eff
        real(rk), intent(in) :: dz_sed

        ! Outputs
        real(rk), intent(out) :: swi_flux       

        ! Locals
        real(rk) :: L_wat, L_sed, logterm
        real(rk) :: Pr, Sc, beta, r_c, G_wat, G_sed, c_swi
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
! Using constant viscosity for now, later compute it from dynamic viscosity. 
        Sc = kin_visc / max(D_sol, eps)

        !--------------------------
        ! Beta function
        ! Approximation for Smooth-bottom only
        ! Appendix B in Umlauf et al. (2023)
        !--------------------------
        beta = 14.8_rk * Sc**(2.0_rk/3.0_rk)

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

        real(rk) :: Fdep

        ! Only downward movement deposits into sediment.
        if (vel_swi < 0._rk) then
            ! Flux
            Fdep = (-vel_swi) * Cw_bot              ! [m/s]*[mass/m3]=[mass/m2/s]
            !Fdep = min(Fdep, Cw_bot * dz_w_bot / dt)   ! prevents immediate negative 
            ! Update concentration at the bottom layer of the water column (Forward Euler)
            Cw_bot = Cw_bot - dt * Fdep / dz_w_bot
            ! update sediment surface concentration (Forward Euler)
            Cbulk_sed_top = Cbulk_sed_top + dt * Fdep / dz_sed_top
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
    !          = alpha(k) * (C_bw - C_s(k)) * Hpw(k) * dt   [mol m^-2]
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
    subroutine apply_bioirrigation(dt, nsed, ntotal, alpha, porewat_thickness, dz_wat_btm, k_wat_btm, concentration)
        real(rk), intent(in)    :: dt
        integer,  intent(in)    :: nsed, ntotal
        real(rk), intent(in)    :: alpha(nsed)               ! Bioirrigation alpha [s^-1]       
        real(rk), intent(in)    :: porewat_thickness(nsed)   ! [m] porewater capacity per m2
        real(rk), intent(in)    :: dz_wat_btm                ! [m] Thickness of the bottom layer in water column        
        integer,  intent(in)    :: k_wat_btm                 ! Indices for water bottom and sediment 
        real(rk), intent(inout) :: concentration(ntotal)     ! Concentration per porewater volume in sediments and water volume in the water column

        integer  :: k
        real(rk) :: Hpw, alpha_k
        real(rk) :: Cbw, Cs, dM, dMtot
        real(rk) :: rmax, lambda
        real(rk), parameter :: eps = 1.0e-30_rk  

        rmax    = maxval(alpha) * dt
        lambda  = dt * sum(alpha * porewat_thickness) / max(dz_wat_btm, eps)
        ! --- Stability / positivity checks for explicit Euler 
        if (rmax   > 1._rk)  stop 'bioirrigation: alpha*dt > 1 (needs more global substeps)'
        if (lambda > 1._rk)  stop 'bioirrigation: lambda > 1 (bottom water can go negative)'

        Cbw = concentration(k_wat_btm)
        dMtot = 0._rk

        do k = 1, nsed
            alpha_k = alpha(k)
            if (alpha_k <= 0._rk) cycle

            ! The amount of mass that can be exchanged depends on how much porewater exists in that layer.
            ! Hpw is the actual porewater volume available to be exchanged
            Hpw = max(porewat_thickness(k), eps)       ! [m] porewater capacity per m2
            Cs = concentration(k)                      ! porewater concentration [mol m-3_pw]

            ! Compute how much mass per m2 is exchanged over dt for each sediment layer
            ! The exchanged volume per unit area is ΔVk = alpha_k* Hpw * dt 
            !   ΔM_k = (Cbw - Ck) * ΔVk
            ! Exchange mass per unit area [mol m-2] over dt
            dM = alpha_k * (Cbw - Cs) * Hpw * dt

            ! Update sediment porewater concentration in the kth sediment layer
            ! by distributing the mass over its porewater volume per m2
            concentration(k) = Cs + dM / Hpw

            ! Accumulate mass removed (or added) from bottom water (opposite sign)
            dMtot = dMtot - dM
        end do

        ! Update bottom water concentration
        concentration(k_wat_btm) = Cbw + dMtot / dz_wat_btm

    end subroutine apply_bioirrigation

end module sediment_exchange