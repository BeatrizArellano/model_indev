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

module turbulence
  use precision_types, only: rk
  use physics_params,  only: rho0, gravity, kappa, mol_nu, PhysicsParams
  use tridiagonal,     only: TridiagCoeff, init_tridiag, solve_tridiag
  implicit none
  private
  public :: TURBULENCE_ke, tke_min

  real(rk), parameter :: tke_min = 3.0e-6_rk   ! Minimum turbulent kinetic energy [m2/s2]
  real(rk), parameter :: eps_min = 5.0e-10_rk  ! Minimum dissipation rate for turbulent kinetic energy
  real(rk), parameter :: L_min   = 0.01_rk
  real(rk), parameter :: cm0     = 0.527_rk    ! neutral-limit constant of the Canuto k–ε closure scheme


contains

    !---------------------------------------------------------------------------
    ! Canuto k-epsilon closure turbulence scheme: builds SS, NN, P, B; solves k and ε;
    ! updates cmue1, cmue2, and finally Nz, Kz. Uses layer interface values at 0..N.
    subroutine TURBULENCE_ke(N, dt, params, h, density, velx, vely, &
                             u_taus, u_taub, z0s, z0b,           &
                             Kz, Nz, tke, eps, Lscale, NN, SS, Ri,          &
                             cmue1, trid, is_first_step)

        integer,               intent(in)    :: N
        real(rk),              intent(in)    :: dt
        type(PhysicsParams),   intent(in)    :: params
        real(rk),              intent(in)    :: h(1:N)
        real(rk),              intent(in)    :: density(1:N)
        real(rk),              intent(in)    :: velx(1:N), vely(1:N)
        real(rk),              intent(in)    :: u_taus, u_taub, z0s, z0b
        real(rk),              intent(inout) :: Kz(0:N), Nz(0:N), Lscale(0:N)  
        real(rk),              intent(inout) :: tke(0:N), eps(0:N), cmue1(0:N)
        real(rk),              intent(inout) :: NN(0:N), SS(0:N), Ri(0:N)          
        type(TridiagCoeff),    intent(inout) :: trid
        logical,               intent(in)    :: is_first_step

        ! Local variables
        real(rk) :: P(0:N), B(0:N), tkeold(0:N)
        real(rk) :: as(0:N), an(0:N)
        real(rk) :: sm(0:N), sh(0:N), cmue2(0:N)
        real(rk) :: x, LLk
        integer  :: i
        real(rk) :: du, dv, dz_imh

        ! Shear (SS), Buoyancy freq (NN), Production P, Buoyancy term B 
        SS(0:N) = 1.0e-6_rk
        NN(0:N) = 0.0_rk
        P (0:N) = 0.0_rk
        B (0:N) = 0.0_rk

        do i=1, N-1            
            ! centered shear from u,v
            dz_imh = 0.5_rk*(h(i)+h(i+1))       ! Distance between centres of layers
            du = velx(i+1) - velx(i)
            dv = vely(i+1) - vely(i)
            ! NOTE: Previously S2P3 had a super long line to compute SS using cnpar but the terms cancelled out, so it is reduced to this
            !SS(i) = SS(i) + (du**2 + dv**2) / (h(i) * h(i+1))
            SS(i) = (du/dz_imh)**2 + (dv/dz_imh)**2

            ! buoyancy frequency N^2 = -(g/rho0) * drho/dz
            NN(i) = -(gravity/rho0) * (density(i+1) - density(i)) / dz_imh

            ! provisional production/buoyancy terms use previous Nz/Kz; 
            ! Provide reasonable initial Nz/Kz
            P(i) = Nz(i) * SS(i)
            B(i) = -Kz(i) * NN(i)
            Ri(i) = NN(i) / max(SS(i), 1.0e-10_rk)  ! Richardson number
        end do

        ! boundaries
        SS(0) = SS(1);   NN(0) = NN(1);   P(0) = P(1);   B(0) = B(1)
        SS(N) = SS(N-1); NN(N) = 0.0_rk; Ri(N) = 0.0_rk; Ri(0) = NN(0)/max(SS(0), 1.0e-10_rk)        

        tkeold = tke  

        ! ---- Solve TKE (k) ----
        call tke_calc(N, dt, h, Nz, P, B, tke, eps, u_taus, u_taub, trid)        

        ! if first time: neutral init to generate an initial cmue1 to pass to dissipation
        
        if (is_first_step) then
            an = 0.0_rk; as = 0.0_rk
            call stability_funcs(N, an, as, cmue1, cmue2, sm, sh)
        end if

        ! ---- Solve epsilon and get length scale ----
        call dissipation(N, dt, h, Nz, tke, tkeold, NN, P, B, z0b, z0s, u_taus, u_taub, eps, cmue1, Lscale, trid)

        ! ---- Stability functions (Canuto): cmue1, cmue2 ----
        do i=0,N
            LLk  = Lscale(i)*Lscale(i) / tke(i)
            as(i)= LLk * SS(i)
            an(i)= LLk * NN(i)
        end do
        
        call stability_funcs(N, an, as, cmue1, cmue2, sm, sh)

        ! ---- Eddy coefficients Nz, Kz  --------
        do i=0,N
            x    = sqrt(tke(i)) * Lscale(i)
            Nz(i)= cmue1(i) * x
            Kz(i)= cmue2(i) * x
        end do

        ! Apply effect of friction and roughness. 
        ! Removed as it is already handled properly in tke_calc and dissipation
        !Nz(0)=kappa*u_taub*z0b; Nz(N)=kappa*u_taus*z0s
        !Kz(0)=kappa*u_taub*z0b; Kz(N)=kappa*u_taus*z0s

        ! Background + limits (no molecular viscosity here)
        do i=0,N
            if (Nz(i) < params%Nz_bg) Nz(i) = params%Nz_bg
            if (Kz(i) < params%Kz_bg) Kz(i) = params%Kz_bg
            !Nz(i) = Nz(i) + mol_nu
            !Kz(i) = Kz(i) + mol_nu
            if (Nz(i) > params%vismax - mol_nu) Nz(i) = params%vismax - mol_nu
            if (Kz(i) > params%vismax - mol_nu) Kz(i) = params%vismax - mol_nu
        end do

    end subroutine TURBULENCE_ke

    !---------------------------------------------------------------------------
    ! Solve TKE (k): tridiagonal implicit step with production/buoyancy split.
    subroutine tke_calc(N, dt, h, Nz, P, B, tke, eps, u_taus, u_taub, tri)
      integer,            intent(in)    :: N
      real(rk),           intent(in)    :: dt
      real(rk),           intent(in)    :: h(1:N), Nz(0:N), P(0:N), B(0:N)
      real(rk),           intent(inout) :: tke(0:N)
      real(rk),           intent(in)    :: eps(0:N), u_taus, u_taub
      type(TridiagCoeff), intent(inout) :: tri

      real(rk) :: avh(0:N)
      real(rk) :: pminus(0:N), pplus(0:N)
      real(rk) :: prod, buoyan, diss, cde
      integer  :: i
    
      cde = cm0**3

      ! interface diffusivity for k
      do i=1,N-1
        avh(i) = 0.5_rk*(Nz(i-1) + Nz(i))
      end do
      avh(1)=0.0_rk; avh(N)=0.0_rk


      ! split sources 
      do i=N-1,1,-1
        prod   = P(i)
        buoyan = B(i)
        diss   = eps(i)
        if (prod+buoyan > 0.0_rk) then
           pplus(i)  =  prod + buoyan
           pminus(i) =  diss
        else
           pplus(i)  =  prod
           pminus(i) =  diss - buoyan
        end if
      end do

      call init_tridiag(tri, N)

      do i=1,N-1
        tri%au(i) = -2.0_rk*dt*avh(i)  /(h(i)+h(i+1))/h(i)
        tri%cu(i) = -2.0_rk*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
        tri%bu(i) =  1.0_rk - tri%au(i) - tri%cu(i) + pminus(i)*dt/tke(i)
        tri%du(i) = (1.0_rk + pplus(i)*dt/tke(i)) * tke(i)
      end do

      ! Solve for i=1..N-1
      call solve_tridiag(1, N-1, tri, tke(1:N-1))

      ! Boundary values 
      tke(0) = u_taub*u_taub/sqrt(cm0*cde)
      tke(N) = u_taus*u_taus/sqrt(cm0*cde)

      ! Minimum value for tke
      where(tke.lt.tke_min) tke=tke_min
    end subroutine tke_calc

    !---------------------------------------------------------------------------
    ! Solve epsilon (ε) and compute length scale L = cde*k^(3/2)/ε.
    subroutine dissipation(N, dt, h, Nz, tke, tkeold, NN, P, B, z0b, z0s, u_taus, u_taub, eps, cmue1, Lscale, tri)
      integer,            intent(in)    :: N
      real(rk),           intent(in)    :: dt
      real(rk),           intent(in)    :: h(1:N), Nz(0:N), tke(0:N), tkeold(0:N), P(0:N), B(0:N)
      real(rk),           intent(in)    :: z0b, z0s, u_taus, u_taub
      real(rk),           intent(inout) :: eps(0:N), cmue1(0:N), NN(0:N)
      real(rk),           intent(out)   :: Lscale(0:N)
      type(TridiagCoeff), intent(inout) :: tri

      real(rk) :: ce1, ce2, galp, cde
      real(rk) :: sig_e(0:N), avh(0:N), flux(0:N)
      real(rk) :: pminus(0:N), pplus(0:N)
      real(rk) :: prod, buoyan, diss, cee3, epslim
      integer :: i

      ce1  = 1.44_rk
      ce2  = 1.92_rk
      galp = 0.27_rk !0.53_rk
      cde  = cm0**3

      !craig_m = sqrt(1.5_rk*cm_craig*cm_craig/(kappa*kappa))
      !sig_e0  = ((4.0_rk/3.0_rk)*craig_m + 1.0_rk ) * (craig_m + 1.0_rk) * (kappa*kappa)/(ce2*cm_craig*cm_craig)
      !sig_e1  =  (kappa*kappa*cm0)/(ce2-ce1)/cde  

      ! Without surface wave contribution
      sig_e(0:N)= (kappa*kappa*cm0)/(ce2-ce1)/cde           !sig_e1 previously


      ! harmonic-like average for eps diffusion coeff
      do i=1,N
         avh(i) = 0.5_rk*(Nz(i-1)/sig_e(i-1) + Nz(i)/sig_e(i))
      end do

      ! Flux boundary contributions (S2P3)
      flux(1:N-1) = 0.0_rk
      flux(1)     = avh(1) * cde * (tkeold(1)**1.5_rk)/ (kappa*(z0b + h(1))**2)      
      ! Without surface wave impact on flux
      flux(N-1)   = cmue1(N-1)*sqrt(tkeold(N-1))*kappa*(h(N)+z0s)/sig_e(N-1)*cde*tkeold(N-1)**1.5/(kappa*(z0s+h(N))**2)    
                                                                        
      avh(1)=0.0; avh(N)=0.0


      ! split sources
      do i=1,N-1
        if (B(i) > 0.0_rk) then
          cee3 = 1.0_rk
        else
          cee3 = -0.629_rk
        end if
        prod   = ce1*eps(i)/tkeold(i) * P(i)
        buoyan = cee3*eps(i)/tkeold(i) * B(i)
        diss   = ce2*eps(i)*eps(i)/tkeold(i)
        if (prod+buoyan > 0.0_rk) then
          pplus(i)  = prod + buoyan
          pminus(i) = diss
        else
          pplus(i)  = prod
          pminus(i) = diss - buoyan
        end if
      end do

      call init_tridiag(tri, N)

      do i=1,N-1
         tri%au(i) = -2.0_rk*dt*avh(i)/(h(i)+h(i+1))/h(i)
         tri%cu(i) = -2.0_rk*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
         tri%bu(i) =  1.0_rk - tri%au(i) - tri%cu(i) + pminus(i)*dt/eps(i)
         tri%du(i) = (1.0_rk + pplus(i)*dt/eps(i))*eps(i) + &
                     flux(i)*dt/(0.5_rk*(h(i)+h(i+1)))
      end do

      call solve_tridiag(1, N-1, tri, eps(1:N-1))         ! Passing only the slice we're interested in solving as boundaries are managed later  

      ! Boundaries
      eps(0) = cde * sqrt(tke(0)*tke(0)*tke(0))/kappa/z0b
      eps(N) = cde * sqrt(tke(N)*tke(N)*tke(N))/kappa/z0s

      ! L scale with stratification cap and floors
      do i=0,N
         if (NN(i) > 0.0_rk) then
            epslim = cde/sqrt(2.0_rk)/galp * tke(i)*sqrt(NN(i))
         else
            epslim = eps_min
         end if
         if (eps(i) < epslim) eps(i) = epslim
         Lscale(i) = cde * sqrt(tke(i)*tke(i)*tke(i))/eps(i)
!
         !if (NN(i) .eq. 0.0_rk) NN = 1.0e-6_rk !NN(i) = 1.0e-6_rk <- Error?
         !if (Lscale(i).gt.(0.267_rk*sqrt(2.0_rk*tke(i)/NN(i)))) then
         !   Lscale(i) =  0.267_rk*sqrt(2.0_rk*tke(i)/NN(i)) 
         !end if
         if (NN(i) > 0.0_rk) then
          ! Only for stable stratification apply an upper cap
          if (Lscale(i) > 0.267_rk*sqrt(2.0_rk*tke(i)/max(NN(i), 1.0e-6_rk))) then
              Lscale(i) = 0.267_rk*sqrt(2.0_rk*tke(i)/max(NN(i), 1.0e-6_rk))
          end if
         end if
         if (Lscale(i) < L_min) Lscale(i) = L_min
      end do
    end subroutine dissipation

    !---------------------------------------------------------------------------
    ! Canuto stability functions: sm, sh and cmue1, cmue2 
    subroutine stability_funcs(N, an, as, cmue1, cmue2, sm, sh)
      integer,  intent(in)  :: N
      real(rk), intent(in)  :: an(0:N), as(0:N)
      real(rk), intent(out) :: cmue1(0:N), cmue2(0:N), sm(0:N), sh(0:N)

      real(rk), parameter :: tnmin = -12.27_rk
      real(rk), parameter :: a2_cm03 = 2.0_rk / cm0**3
                               
      real(rk), parameter :: d(6) = [4.2483e+2_rk, 2.7132e+1_rk, 3.0498889_rk, 2.304e-1_rk, 1.3866e-1_rk, -8.953856e-4_rk]
      real(rk), parameter :: s(6) = [2.2728405e+1_rk, 9.245632e-1_rk, -6.42e-3_rk, 2.38e+1_rk, 2.4e-1_rk, 4.6702411e-2_rk]     


      real(rk) :: tn, ts, tsmax, dd
      integer  :: i

      do i=1, N-1
         tn = 4.0_rk / cm0**6 * an(i) 
         if (tn < tnmin) tn = tnmin
         ts = 4.0_rk / cm0**6 * as(i)
         tsmax = (d(1)+d(2)*tn + d(4)*tn*tn) / (d(3) + d(5)*tn)
         if (ts > tsmax) ts = tsmax

         dd   = d(1)+d(2)*tn + d(3)*ts + d(4)*tn*tn + d(5)*tn*ts + d(6)*ts*ts
         sm(i)= (s(1)+s(2)*tn + s(3)*ts) / dd
         sh(i)= (s(4)+s(5)*tn + s(6)*ts) / dd
         cmue1(i) = a2_cm03 * sm(i)
         cmue2(i) = a2_cm03 * sh(i)
      end do

      ! Boundaries
      cmue1(0)=cmue1(1); cmue1(N)=cmue1(N-1)
      cmue2(0)=cmue2(1); cmue2(N)=cmue2(N-1)
      !sm(0)=sm(1); sm(N)=sm(N-1)
      !sh(0)=sh(1); sh(N)=sh(N-1)
    end subroutine stability_funcs

end module turbulence
