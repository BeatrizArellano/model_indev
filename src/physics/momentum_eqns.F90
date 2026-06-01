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
    use physics_params,  only: rho_air, rho0, kappa, mol_nu, gravity, z0s_min
    use precision_types, only: rk    
    use trigonometrics,  only: deg2rad
    implicit none
    private

    public :: EQN_PRESSURE, EQN_CORIOLIS
    public :: compute_surface_stress, compute_bottom_stress

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

    pure subroutine compute_surface_stress(u10, v10, charnock, tau_x, tau_y, u_taus, z0s, windspeed)
                                            
        implicit none
        real(rk), intent(in)     :: u10, v10           ! 10m wind components [m s-1]
        real(rk), intent(in)     :: charnock
        real(rk), intent(out)    :: tau_x, tau_y       ! surface stress [N m-2]
        real(rk), intent(out)    :: u_taus, z0s
        real(rk), intent(out)    :: windspeed          ! optional: return |U10|     
        

        ! Locals
        real(rk) :: tau_mag, Cd, fac  
        
        ! Magnitude of 10m windspeed
        windspeed = sqrt(u10*u10 + v10*v10)

        ! Defaults (important for intent(out))
        tau_x  = 0._rk
        tau_y  = 0._rk
        u_taus = 0._rk
        z0s    = z0s_min

        ! Calm case
        if (windspeed <= 1.0e-10_rk) return

        ! Smith & Banke (1975): Cd = (0.63 + 0.066*wind_speed)×1e-3 (In Simon and Sharples)
        Cd  = (0.63_rk + 0.066_rk*windspeed) * 1e-3_rk

        ! Air-side stress exerted on ocean:
        ! Stress τ = rho_air * Cd * |ws| * u (or v)
        fac    = rho_air * Cd * windspeed
        tau_x = fac * u10
        tau_y = fac * v10

        !	surface stress boundary conditions...   
        tau_mag = sqrt(tau_x*tau_x + tau_y*tau_y)
        !u_taus: water-side friction velocity or shear velocity
        u_taus  = sqrt(tau_mag / rho0)
        ! water-side roughness length -> Formulations from NEMO have found a charnock value of 1400 is appropriate (e.g. Alari et al., 2016)
        z0s     = max(charnock * (u_taus*u_taus) / gravity, z0s_min)    
    end subroutine compute_surface_stress


    subroutine compute_bottom_stress(u, v, dz_btm, h0b, tau_bx, tau_by, u_taub, z0b, stressb)
        ! Inputs
        real(rk), intent(in)  :: u(:), v(:)      ! velocity profiles (bottom-to-top)
        real(rk), intent(in)  :: dz_btm          ! bottom layer thickness [m]
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
        z0r = 0.033_rk * h0b
        z0b = max(z0r, 1.0e-9_rk)

        ! Speed magnitude at bottom cell centre
        Uc = sqrt(u(1)*u(1) + v(1)*v(1))

        if (Uc <= 1.0e-10_rk) then
            ! essentially no flow so no quadratic stress
            return
        end if

        ! ---- Bottom friction: using law-of-the-wall to compute friction velocity
        ! iterate a few times because z0sm (smooth term) depends on ustar via rb*Uc, both terms depend on each other
        do it=1,3
            ! log-law factor at the actual sampling height zc = h(1)/2 and keeping it out of the viscous sublayer (10*z0b)
            logarg = log((z0b + zc) / z0b)                     ! Log part
            if (logarg <= 1.0e-12_rk) error stop "compute_bottom_stress: invalid logarg"
            rb     = kappa / logarg                            ! rb = ustar/U(z) = kappa/ln(z/z0)                    
            u_taub= rb * Uc                                    ! shear velocity from log-law
            if (mol_nu>0.0_rk .and. u_taub>0.0_rk) then
                z0sm = 0.11_rk * mol_nu/max(1.0e-15_rk, u_taub)
            else
                z0sm = 0.0_rk
            end if
            ! updated roughness length 
            ! Colebrook-White style transition formula 
            ! Include smooth and rough limits
            z0new = max(z0sm + z0r, 1.0e-15_rk)    
            z0b = z0new
        end do

        ! stress magnitude
        stressb = rho0 * u_taub * u_taub  ! Using rho0 to keep a Boussinesq-style

        ! direction: oppose the flow
        tau_bx = -stressb * (u(1) / Uc)
        tau_by = -stressb * (v(1) / Uc)

    end subroutine compute_bottom_stress

    !----------------------------------------------------------------------------------
    ! -------- Legacy subroutines
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
    
end module momentum_eqns