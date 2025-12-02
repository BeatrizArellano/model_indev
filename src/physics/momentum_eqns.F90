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

    public :: EQN_PRESSURE, EQN_CORIOLIS, wind_stress_from_uv, EQN_FRICTION, update_surface_friction

contains
    ! Equation of motion - pressure term
    ! Uniform barotropic pressure-gradient acceleration applied to all levels
    pure subroutine EQN_PRESSURE(dt, Pgrad, u)
        real(rk), intent(in)    :: dt          ! time step [s]
        real(rk), intent(in)    :: Pgrad       ! acceleration [m s-2] (x or y component)
        real(rk), intent(inout) :: u(:)        ! velocity profile for the given component [m s-1]
        ! Equivalent to S2P3 loop: u(i) = u(i) + dt * Pgrad for i=1..N
        u = u + dt * Pgrad
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



    !**********************************************************
    !       Equation of motion subroutine - friction term
    !**********************************************************
    pure subroutine EQN_FRICTION(vel_comp_old, vel_comp_new, vel_comp2_bottom, Nz, h, dt, h0b, density, &
                                 tau_surf, u_taub, z0b, stressb)
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
        real(rk), intent(in)  :: vel_comp_old(:), vel_comp2_bottom, Nz(:), h(:), density(:), dt, h0b, tau_surf
        ! Outputs
        real(rk), intent(out) :: vel_comp_new(size(vel_comp_old)), u_taub, z0b, stressb

        integer  :: i, N
        real(rk) :: dz_imh, dz_iph, flux_dn, flux_up, dv
        real(rk) :: zc, Uc, rb, ustar_b, z0r, z0sm, z0new, tau_mag, logarg
        real(rk) :: rho_bed
        integer  :: it

        N     = size(vel_comp_old); vel_comp_new = vel_comp_old
        zc    = 0.5_rk*h(1)

        rho_bed  = density(1)

        ! ---- Interior diffusion of horizontal momentum
        do i = 2, N-1
            ! distances between centres across layer interfaces
            dz_imh = 0.5_rk*(h(i-1) + h(i))      ! face i-1/2
            dz_iph = 0.5_rk*(h(i)   + h(i+1))    ! face i+1/2

            ! fluxes (flux_up - flux_dn)
            ! flux leaving cell i downward into cell i−1
            flux_dn = Nz(i-1) * (vel_comp_old(i)   - vel_comp_old(i-1)) / dz_imh   ! at i-1/2
            ! flux leaving cell i upward into cell i+1
            flux_up = Nz(i  ) * (vel_comp_old(i+1) - vel_comp_old(i  )) / dz_iph   ! at i+1/2

            dv = dt * ((flux_up - flux_dn) / h(i))   ! tendency
            vel_comp_new(i) = vel_comp_old(i) + dv
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
        z0r = 0.03_rk*h0b           ! converts Nikuradse height (h0b) to roughness length (approx k_s/33).
        z0b = max(z0r, 1.0e-6_rk)   ! initial roughness (rough limit only to start)

        ! iterate a few times because z0sm (smooth term) depends on ustar via rb*Uc, both terms depend on each other
        do it=1,3
            ! log-law factor at the actual sampling height zc = h(1)/2 and keeping it out of the viscous sublayer (10*z0b)
            logarg = log(max(zc, 10.0_rk*z0b) / z0b)           ! Log part
            rb     = kappa / logarg                            ! rb = ustar/U(z) = kappa/ln(z/z0)
            ! Magnitude of speed at the centre of bottom layer
            Uc     = sqrt(vel_comp_old(1)*vel_comp_old(1) + vel_comp2_bottom*vel_comp2_bottom)            
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
        flux_up = 0.0_rk
        if (N >= 2) then
            dz_iph  = 0.5_rk*(h(1) + h(2))                                
            flux_up = Nz(1) * (vel_comp_old(2) - vel_comp_old(1)) / dz_iph          ! upward diffusive flux
        end if
        ! quadratic bottom sink with height-consistent coefficient rb^2
        vel_comp_new(1) = vel_comp_old(1) + dt * (flux_up / h(1) - (rb*rb/h(1))*Uc*vel_comp_old(1))

        ! Bottom stress boundary conditiosn
        stressb    = rho_bed * ustar_b * ustar_b
        u_taub     = ustar_b

        ! ---- Surface: wind stress
        flux_dn = 0.0_rk
        if (N >= 2) then
            dz_imh  = 0.5_rk*(h(N-1) + h(N))                 ! face N-1/2 spacing
            flux_dn = Nz(N-1) * (vel_comp_old(N) - vel_comp_old(N-1)) / dz_imh
        end if
        dv = dt * (-flux_dn / h(N) + tau_surf/(rho0*h(N)))  ! -> rho0 or density at the surface?
        vel_comp_new(N) = vel_comp_old(N) + dv

    end subroutine EQN_FRICTION
    
end module momentum_eqns