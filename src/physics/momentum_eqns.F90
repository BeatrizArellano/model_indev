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

module momentum_eqns
    use precision_types, only: rk
    use physics_params,  only: rho_air, rho0, kappa, mol_nu, gravity, z0s_min
    use trigonometrics,  only: deg2rad
    implicit none
    private

    public :: EQN_PRESSURE, EQN_CORIOLIS, wind_stress_from_uv, update_surface_friction
    public :: compute_bottom_stress

contains

    ! Equation of motion - pressure term
    ! Uniform barotropic pressure-gradient acceleration applied to all levels
    pure subroutine EQN_PRESSURE(dt, Pgrad_x, Pgrad_y, u, v)
        real(rk), intent(in)    :: dt                   ! time step [s]
        real(rk), intent(in)    :: Pgrad_x, Pgrad_y     ! acceleration [m s-2] (x or y component)
        real(rk), intent(inout) :: u(:), v(:)           ! velocity profile for the given component [m s-1]


        u = u + dt * Pgrad_x
        v = v + dt * Pgrad_y
    end subroutine EQN_PRESSURE

    !**********************************************************
    !       Semi-implicit Coriolis.
    !**********************************************************
    pure subroutine EQN_CORIOLIS(u, v, timestep, f0)
        real(rk), intent(inout) :: u(:), v(:)
        real(rk), intent(in)    :: timestep, f0
        real(rk), allocatable   :: u0(:), v0(:)
        real(rk) :: alf, alf1, alf2
        
        alf  = f0 * timestep                      ! constants used for semi-implicit Coriolis
        alf1 = 1.0_rk / (1.0_rk + 0.25_rk*alf*alf)
        alf2 = 1.0_rk - (0.25_rk*alf*alf)             

        u0 = u;  v0 = v                            ! keeps the old state for updating both u and v
        u  = alf1 * ( alf2*u0 + alf*v0)
        v  = alf1 * (-alf*u0 + alf2*v0)
    end subroutine EQN_CORIOLIS


    pure subroutine wind_stress_from_speed_dir(windspeed, wind_direction, stressx, stressy)
        implicit none
        real(rk), intent(in)  :: windspeed, wind_direction
        real(rk), intent(out) :: stressx, stressy      ! surface stress [N m^-2]
        real(rk) :: th, Cd, u10, v10, fac

        th   =  deg2rad(wind_direction)
        u10  = -windspeed * sin(th)      ! East component
        v10  = -windspeed * cos(th)      ! north
        ! Smith & Banke (1975): Cd = (0.63 + 0.066*wind_speed)×1e-3 (In Simon and Sharples)
        Cd = (0.63_rk + 0.066_rk*windspeed) * 1e-3_rk 
          ! Stress τ = rho_air * Cd * |ws| * u (or v)
        fac   = rho_air * Cd * windspeed
        stressx = fac * u10
        stressy = fac * v10
        ! w_stress=sqrt(stressx*stressx + stressy*stressy) ! also called stresss in S2P3 not needed for the calculations
    end subroutine wind_stress_from_speed_dir

    pure subroutine wind_stress_from_uv(u10, v10, windspeed, stressx, stressy)
        implicit none
        real(rk), intent(in)     :: u10, v10           ! 10m wind components [m s-1]
        real(rk), intent(out)    :: stressx, stressy   ! surface stress [N m-2]
        real(rk), intent(out)    :: windspeed          ! optional: return |U10|
        real(rk) :: Cd, fac

        ! Magnitude of 10m windspeed
        windspeed = sqrt(u10*u10 + v10*v10)

        ! Calm case
        if (windspeed <= 0._rk) then
            stressx = 0._rk
            stressy = 0._rk
            return
        end if
        ! Smith & Banke (1975): Cd = (0.63 + 0.066*wind_speed)×1e-3 (In Simon and Sharples)
        Cd  = (0.63_rk + 0.066_rk*windspeed) * 1e-3_rk
        ! Stress τ = rho_air * Cd * |ws| * u (or v)
        fac    = rho_air * Cd * windspeed
        stressx = fac * u10
        stressy = fac * v10
    end subroutine wind_stress_from_uv

    pure subroutine update_surface_friction(tau_x, tau_y, rho_surf, charnock, u_taus, z0s)
        !  u_taus,           [out] real     friction velocity at the surface u_* [m s-1]
        !  z0s               [out] real     effective surface roughness length  [m]
        implicit none
        real(rk), intent(in)  :: tau_x, tau_y, rho_surf, charnock
        real(rk), intent(out) :: u_taus, z0s
        real(rk) :: tau_mag

        !	surface stress boundary conditions...   
        tau_mag = sqrt(tau_x*tau_x + tau_y*tau_y)
        !u_taus: water-side friction velocity or shear velocity
        u_taus  = sqrt(tau_mag / rho_surf)
        ! water-side roughness length -> Formulations from NEMO have found a charnock value of 1400 is appropriate (e.g. Alari et al., 2016)
        z0s     = max(charnock * (u_taus*u_taus) / gravity, z0s_min)
    end subroutine update_surface_friction



    subroutine compute_bottom_stress(u, v, dz_btm, rho_bed, h0b, tau_bx, tau_by, u_taub, z0b, stressb)
        ! Inputs
        real(rk), intent(in)  :: u(:), v(:)      ! velocity profiles (bottom-to-top)
        real(rk), intent(in)  :: dz_btm          ! bottom layer thickness [m]
        real(rk), intent(in)  :: rho_bed         ! density at bottom cell [kg/m3]
        real(rk), intent(in)  :: h0b             ! Nikuradse roughness height k_s [m]

        ! Outputs
        real(rk), intent(out) :: tau_bx, tau_by  ! bed stress components [N/m2]
        real(rk), intent(out) :: u_taub          ! bed friction velocity [m/s]
        real(rk), intent(out) :: z0b             ! effective roughness length [m]
        real(rk), intent(out) :: stressb         ! bed stress magnitude [N/m2]

        ! Locals
        real(rk) :: zc, Uc, logarg
        real(rk) :: z0r, rb, z0sm, z0new
        integer  :: it

        ! Default outputs
        tau_bx  = 0.0_rk
        tau_by  = 0.0_rk
        u_taub  = 0.0_rk
        stressb = 0.0_rk

        ! Sampling height at centre of bottom layer
        zc = 0.5_rk * dz_btm

        ! Convert Nikuradse height (h0b) to roughness length (approx k_s/33).
        z0r = 0.03_rk * h0b
        z0b = max(z0r, 1.0e-9_rk)

        ! Speed magnitude at bottom cell centre
        Uc = sqrt(u(1)*u(1) + v(1)*v(1))

        if (Uc <= 1.0e-10_rk) then
            ! essentially no flow so no quadratic stress
            return
        end if

        ! iterate a few times because z0sm (smooth term) depends on ustar via rb*Uc, both terms depend on each other
        do it=1,3
            ! log-law factor at the actual sampling height zc = h(1)/2 and keeping it out of the viscous sublayer (10*z0b)
            logarg = log((z0b + zc) / z0b)                     ! Log part
            if (logarg <= 1.0e-12_rk) error stop "compute_bottom_stress: invalid logarg"
            rb     = kappa / logarg                            ! rb = ustar/U(z) = kappa/ln(z/z0)                    
            u_taub= rb * Uc                                    ! shear velocity from log-law
            if (mol_nu>0.0_rk .and. u_taub>0.0_rk) then
                z0sm = 0.1_rk * mol_nu/max(mol_nu, u_taub)
            else
                z0sm = 0.0_rk
            end if
            ! updated roughness length 
            z0new = max(z0sm + z0r, 1.0e-9_rk)    
            z0b = z0new
        end do

        ! stress magnitude
        stressb = rho_bed * u_taub * u_taub

        ! direction: oppose the flow
        tau_bx = -stressb * (u(1) / Uc)
        tau_by = -stressb * (v(1) / Uc)

    end subroutine compute_bottom_stress


    !**********************************************************
    !       Equation of motion subroutine - friction term
    !**********************************************************
    subroutine EQN_FRICTION(u_old, v_old, u_new, v_new, Nz, h, dt, h0b, density,  &
                                 tau_x, tau_y, u_taub, z0b, stressb)
        !  vel_comp_old,     [in ] real(:)  bottom-to-top component to update (u or v) at layers' centres [m s-1]
        !  vel_comp_new,     [out] real(:)  updated component after diffusion + bottom drag + surface stress [m s-1]
        !  vel_comp2_bottom  [in ] real     other component at i=1 (bottom-layer centre) used to form |U| [m s-1]
        !  Nz,               [in ] real(:)  eddy viscosity (momentum) at interfaces (0..N) [m^2 s-1]
        !  h,                [in ] real(:)  layer thicknesses (centres 1..N) [m]
        !  dt,               [in ] real     time step [s]
        !  h0b,              [in ] real     Nikuradse bed roughness height k_s (physical roughness) [m]
        !  density                 real     Density 
        !  tau_surf,         [in ] real     surface stress component (τ_sx or τ_sy) applied at top layer [N m-2]
        !  u_taub,           [out] real     friction velocity at bed u_* [m s-1]        
        !  z0b              [out] real     effective bottom roughness length  [m]

        implicit none
        ! Inputs
        real(rk), intent(in)  :: u_old(:), v_old(:)
        real(rk), intent(in)  :: Nz(:), h(:), density(:)
        real(rk), intent(in)  :: dt, h0b, tau_x, tau_y
        real(rk), intent(out) :: u_new(size(u_old)), v_new(size(v_old))
        real(rk), intent(out) :: u_taub, z0b, stressb

        integer  :: i, N, it
        real(rk) :: dz_imh, dz_iph, flux_dn_u, flux_up_u, flux_dn_v, flux_up_v
        real(rk) :: dv_u, dv_v
        real(rk) :: zc, Uc, rb, ustar_b, z0r, z0sm, z0new, logarg
        real(rk) :: rho_bed, rho_surf

        N = size(u_old)
        u_new = u_old
        v_new = v_old

        rho_bed  = density(1)
        rho_surf = density(N)

        ! ---- Interior diffusion of horizontal momentum
        do i = 2, N-1
            ! distances between centres across layer interfaces
            dz_imh = 0.5_rk*(h(i-1) + h(i))      ! face i-1/2
            dz_iph = 0.5_rk*(h(i)   + h(i+1))    ! face i+1/2

            ! fluxes (flux_up - flux_dn)
            ! flux leaving cell i downward into cell i−1
            flux_dn_u = Nz(i-1) * (u_old(i)   - u_old(i-1)) / dz_imh   ! at i-1/2
            flux_dn_v = Nz(i-1) * (v_old(i)   - v_old(i-1)) / dz_imh
            ! flux leaving cell i upward into cell i+1
            flux_up_u = Nz(i  ) * (u_old(i+1) - u_old(i  )) / dz_iph   ! at i+1/2            
            flux_up_v = Nz(i  ) * (v_old(i+1) - v_old(i  )) / dz_iph

            dv_u = dt * ((flux_up_u - flux_dn_u) / h(i))   ! tendency
            dv_v = dt * ((flux_up_v - flux_dn_v) / h(i))
            u_new(i) = u_old(i) + dv_u
            v_new(i) = v_old(i) + dv_v
        end do


        ! ---- Bottom friction: using law-of-the-wall to compute friction velocity
        !   NOTES:
        !   - S2P3 used a constant kb with a 1 m "bed" velocity (via bed_factor).
        !   - That relates stress to a fixed reference height and to a fixed layer thickness.
        !   - Here we compute a height-consistent drag from the log-law at the actual sampling height
        !   Law of the wall: U(z)=(ustar/kappa)*ln(z/z0)
        !   where ustar is friction velocity (shear), z is depth, and z0 is roughness length. 
        !   U(z) is the mean velocity at height z above the bed.
        !   It assumes we are in the logarithmic region of the velocity profile, but at a height larger than z0. 
        !     zc = h1/2. The momentum sink still enters the bottom cell like in S2P3.
        zc  = 0.5_rk*h(1)
        ! Magnitude of speed at the centre of bottom layer
        Uc  = sqrt(u_old(1)*u_old(1) + v_old(1)*v_old(1))

        z0r = 0.03_rk * h0b           ! converts Nikuradse height (h0b) to roughness length (approx k_s/33).
        z0b = max(z0r, 1.0e-6_rk)     ! initial roughness (rough limit only to start)

        ! iterate a few times because z0sm (smooth term) depends on ustar via rb*Uc, both terms depend on each other
       do it=1,3
            ! log-law factor at the actual sampling height zc = h(1)/2 and keeping it out of the viscous sublayer (10*z0b)
            logarg = log(max(zc, 10.0_rk*z0b) / z0b)           ! Log part
            rb     = kappa / logarg                            ! rb = ustar/U(z) = kappa/ln(z/z0)                    
            ustar_b= rb * Uc                                   ! shear velocity from log-law
            if (mol_nu>0.0_rk .and. ustar_b>0.0_rk) then
                z0sm = 0.1_rk*mol_nu/ustar_b
            else
                z0sm = 0.0_rk
            end if
            ! updated roughness length 
            z0new = max(z0sm + z0r, 1.0e-6_rk)    
            z0b = z0new
       end do

        ! diffusive flux from layer 2 across the bottom face (1/2), using centre-to-centre distance
        flux_up_u = 0.0_rk
        flux_up_v = 0.0_rk
        if (N >= 2) then
            dz_iph  = 0.5_rk*(h(1) + h(2))                         
            flux_up_u = Nz(1) * (u_old(2) - u_old(1)) / dz_iph          ! upward diffusive flux
            flux_up_v = Nz(1) * (v_old(2) - v_old(1)) / dz_iph
        end if
        ! quadratic bottom sink with height-consistent coefficient rb^2
        u_new(1) = u_old(1) + dt * (flux_up_u / h(1) - (rb*rb/h(1)) * Uc * u_old(1))
        v_new(1) = v_old(1) + dt * (flux_up_v / h(1) - (rb*rb/h(1)) * Uc * v_old(1))

        ! Bottom stress boundary conditiosn
        stressb    = rho_bed * ustar_b * ustar_b
        u_taub     = ustar_b


        ! ---- Surface: wind stress
        flux_dn_u = 0.0_rk
        flux_dn_v = 0.0_rk
        if (N >= 2) then
            dz_imh  = 0.5_rk*(h(N-1) + h(N))                 ! face N-1/2 spacing
            flux_dn_u = Nz(N-1) * (u_old(N) - u_old(N-1)) / dz_imh
            flux_dn_v = Nz(N-1) * (v_old(N) - v_old(N-1)) / dz_imh
        end if
        ! -> rho0 or density at the surface?
        u_new(N) = u_old(N) + dt * (-flux_dn_u/ h(N) + tau_x/(rho_surf*h(N)))
        v_new(N) = v_old(N) + dt * (-flux_dn_v/ h(N) + tau_y/(rho_surf*h(N)))
    end subroutine EQN_FRICTION
    
end module momentum_eqns