module numerical_stability
  use bio_types,       only: BioEnv
  use precision_types, only: rk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

  implicit none
  private
  public :: compute_phys_subcycles, compute_bio_substeps



contains


  pure subroutine compute_phys_subcycles(h, Nz, Kz, vismax, dt_main, n_sub, dt_sub, ierr, errmsg)
    !   Calculates an integer number of substeps (n_sub) and the substep size (dt_sub)
    !   so explicit vertical diffusion of MOMENTUM is stable over dt_sub.
    !   This is done by calculating the stability condition (CFL)
    !   so that CFL: v*dt_sub/Δz^2 ≤ 1/2 everywhere.
    !
    ! Arrays:
    !   h(1..N)      : layer thickness [m]
    !   Nz(0..N)     : momentum viscosity at interfaces [m2/s]
    !   Kz(0..N)     : scalar diffusivity at interfaces [m2/s]
    !
    !   vismax       : maximum viscosity or diffusivity allowed int he model
    !
    ! Time:
    !   dt_main      : outer step [s]
    !
    ! Outputs:
    !   n_sub      : number of substeps
    !   dt_sub    : inner time-step
    !
    real(rk), intent(in)  :: h(:), Nz(:), Kz(:)
    real(rk), intent(in)  :: dt_main, vismax
    integer,  intent(out) :: n_sub
    real(rk), intent(out) :: dt_sub
    integer,  intent(out) :: ierr
    character(len=*), intent(out) :: errmsg

    real(rk), parameter :: SAFETY      = 0.49_rk      ! inside parabolic limit with momentum diffusion explicit and Nz fixed over the subcycle
    real(rk), parameter :: DT_MIN      = 0.1_rk       ! [s] Stop below this time-step (Physics)
    real(rk), parameter :: D_MIN       = 1.0e-12_rk   ! Minimum value for eddy coefficients so we don't divide by zero
    real(rk), parameter :: BOOST = 15._rk, FRAC_vismax = 0.7_rk

    integer  :: N, i
    real(rk) :: dz, nz_sens, kz_sens
    real(rk) :: dtmin_mom, dtmin_sca, dt_safe

    N = size(h)
    if (size(Nz) /= N+1 .or. size(Kz) /= N+1) then
      ierr = 3; errmsg = 'compute_phys_subcycles: Nz and Kz size must be 0..N (N+1).'
      n_sub = 1; dt_sub = dt_main
      return
    end if


    ! Find explicit diffusion limits on interfaces (1..N-1)
    dtmin_mom = huge(1.0_rk)
    dtmin_sca = huge(1.0_rk)



    do i = 1, N-1
      dz = 0.5_rk * (h(i) + h(i+1))
      nz_sens = max(D_MIN,Nz(i))
      kz_sens = max(D_MIN,Kz(i))

      if (nz_sens <= FRAC_vismax *vismax) then
         nz_sens = min(BOOST * nz_sens, vismax)
      else 
         nz_sens = vismax
      end if
      if (kz_sens <= FRAC_vismax *vismax) then
         kz_sens = min(BOOST*kz_sens, vismax)
      else
        kz_sens = vismax
      end if

      dtmin_mom = min(dtmin_mom, (dz*dz) / (2.0_rk * nz_sens)) 
      dtmin_sca = min(dtmin_sca, (dz*dz) / (2.0_rk * kz_sens))
    end do

    ! Strictest limit across momentum & scalar
    dt_safe = SAFETY * min(dtmin_mom, dtmin_sca)

    ! Integer subcycles so each substep <= dt_safe
    n_sub  = max(1, ceiling(dt_main / dt_safe))
    dt_sub = dt_main / real(n_sub, rk)

    ierr = 0; errmsg = ''
    if (dt_sub < DT_MIN) then
      ierr = 1
      errmsg = 'compute_phys_subcycles: internal time-step < 0.1 s;  Increase layer thickness in the vertical grid.'
    end if
    if (.not. (dt_safe > 0._rk) .or. ieee_is_nan(dt_safe)) then
      ierr = 2
      errmsg = 'compute_phys_subcycles: dt_safe not positive (bad viscosity/diffusivity inputs).'
      n_sub = 1
      dt_sub = dt_main
      return
    end if
  end subroutine compute_phys_subcycles


  ! Computes the number of subcycles required while updating biogeochemistry, considering potential limiting factors for all processes. 
  subroutine compute_bio_substeps(BE, dt_main, frac_max, dt_min, &
                                  nsub_adv, nsub_diff, nsub_bio, &
                                  nsub, dt_sub,                  &
                                  nsub_burial, nsub_bioirr, nsub_dep, nsub_sed_diff)

    !--- Inputs ---
    type(BioEnv), intent(in) :: BE
    real(rk),     intent(in) :: dt_main
    real(rk),     intent(in) :: frac_max
    real(rk),     intent(in) :: dt_min

    !--- Outputs (core) ---
    integer,      intent(out) :: nsub_adv, nsub_diff, nsub_bio
    integer,      intent(out) :: nsub
    real(rk),     intent(out) :: dt_sub

    !--- Optional extra diagnostics (useful for logging) ---
    integer,      intent(out), optional :: nsub_burial, nsub_bioirr, nsub_dep, nsub_sed_diff

    !--- Locals ---
    integer  :: kwb, kws, nwat, nsed
    real(rk) :: cfl_adv, cfl_adv_sol, cfl_adv_sld
    real(rk) :: dt_dummy
    real(rk) :: nsub_real_limit
    integer  :: nsub_limit

    integer :: ns_burial, ns_bioirr, ns_dep, ns_sed_diff, ns_sed_biot, ns_burial_sol, ns_burial_sld


    ! defaults
    nsub_adv  = 1
    nsub_diff = 1
    nsub_bio  = 1
    ns_burial   = 1
    ns_bioirr   = 1
    ns_dep      = 1
    ns_sed_diff = 1
    ns_sed_biot = 1
    ns_burial_sol = 1
    ns_burial_sld = 1
    cfl_adv_sol   = 0._rk
    cfl_adv_sld   = 0._rk


    kwb  = BE%k_wat_btm
    kws  = BE%k_wat_sfc
    nwat = BE%nwat
    nsed = BE%nsed

    !----------------------------------------------------------------
    !  Substeps from vertical transport (sinking or rising)
    !   Only in the Water-column 
    !----------------------------------------------------------------
    call compute_transport_substeps(BE%vel_faces(kwb-1:kws, :), BE%wat_grid%dz(1:nwat), &
                                    dt_main, cfl_adv, nsub_adv)

    !----------------------------------------------------------------
    ! Water-column diffusion stability 
    !----------------------------------------------------------------
    call compute_diffusion_substeps(BE%BS%vert_diff(kwb-1:kws), BE%wat_grid%dz(1:nwat), &
                                    dt_main, BE%params%cnpar, nsub_diff, dt_dummy)

    !----------------------------------------------------------------
    ! Reaction / source integration stability (explicit Euler)
    !----------------------------------------------------------------
    call compute_bioint_substeps(BE, BE%tendency_int, BE%tendency_sf, BE%tendency_bt, &
                                 dt_main, frac_max, nsub_bio, dt_dummy)

    !----------------------------------------------------------------
    ! Sediment-related explicit operators (only if enabled)
    !----------------------------------------------------------------
    if (BE%params%sediments_enabled) then

      ! Sediment burial CFL (solutes + solids velocities).
      call compute_transport_substeps_sed(BE%SED%vel_solutes(0:nsed), BE%sed_grid%dz(1:nsed), &
                                          dt_main, cfl_adv_sol, ns_burial_sol )

      call compute_transport_substeps_sed(BE%SED%vel_solids(0:nsed), BE%sed_grid%dz(1:nsed), &
                                          dt_main, cfl_adv_sld, ns_burial_sld )

      ns_burial = max(ns_burial_sol, ns_burial_sld)


      ! Bioirrigation explicit stability (rmax + lambda)
      call compute_bioirrigation_substeps(alpha             = BE%SED%bioirr(1:nsed), &
                                          porewat_thickness = BE%SED%porewat_thickness(1:nsed), &
                                          dz_wat_btm        = BE%wat_grid%dz(1), &
                                          dt_main           = dt_main, &
                                          nsub              = ns_bioirr)

      ! Particulate deposition explicit positivity
      call compute_deposition_substeps(BE, kwb, BE%wat_grid%dz(1), dt_main, ns_dep)

      ! Sediment diffusion CFL 
      ! If cnpar_sed==1, this is unconditionally stable and can be ignored.
      if (BE%SED%params_SI%cnpar_sed < 1.0_rk) then
        call compute_diffusion_substeps_sed(BE%SED%diff_sed_max(0:nsed), BE%sed_grid%dz(1:nsed), BE%SED%porewat_thickness(1:nsed), &
                                            dt_main, cnpar=BE%SED%params_SI%cnpar_sed,   &
                                            nsub=ns_sed_diff, dt_sub=dt_dummy)

        call compute_diffusion_substeps_sed(BE%SED%Db_eff_solids(0:nsed), BE%sed_grid%dz(1:nsed), BE%SED%solid_thickness(1:nsed), &
                                            dt_main, cnpar=BE%SED%params_SI%cnpar_sed,   &
                                            nsub=ns_sed_biot, dt_sub=dt_dummy)
      else
        ns_sed_diff = 1
        ns_sed_biot = 1
      end if

    end if

    !----------------------------------------------------------------
    ! Global number of substeps
    !----------------------------------------------------------------
    if (BE%params%sediments_enabled) then
      nsub = max(nsub_adv, max(nsub_diff, nsub_bio) )
      nsub = max(nsub, max(ns_burial, max(ns_bioirr, max(ns_dep, max(ns_sed_diff, ns_sed_biot)))) )
    else
      nsub = max(nsub_adv, max(nsub_diff, nsub_bio) )
    end if

    if (nsub < 1) nsub = 1
    dt_sub = dt_main / real(nsub, rk)

    !----------------------------------------------------------------
    ! dt_sub >= dt_min (your existing clamp)
    !----------------------------------------------------------------
    if (dt_min > 0._rk .and. dt_sub < dt_min) then
      nsub_real_limit = dt_main / dt_min
      nsub_limit      = int(nsub_real_limit)  ! floor

      if (nsub_limit < 1) then
        nsub   = 1
        dt_sub = dt_main
      else
        if (nsub > nsub_limit) nsub = nsub_limit
        if (nsub < 1) nsub = 1
        dt_sub = dt_main / real(nsub, rk)
      end if
    end if

    !----------------------------------------------------------------
    ! Optional diagnostics
    !----------------------------------------------------------------
    if (present(nsub_burial))   nsub_burial   = ns_burial
    if (present(nsub_bioirr))   nsub_bioirr   = ns_bioirr
    if (present(nsub_dep))      nsub_dep      = ns_dep
    if (present(nsub_sed_diff)) nsub_sed_diff = max(ns_sed_diff, ns_sed_biot)

  end subroutine compute_bio_substeps



  !--------------------------------------------------------------------
  ! Compute maximum vertical Courant number and numebr of substeps.
  ! This guarantees that each tracer never moves more than one layer each substep
  !
  !  - w_face(0:N) : vertical velocity at interfaces [m s-1]
  !  - dz(1:N)     : layer thickness [m]
  !  - dt          : time step [s]
  !
  !  cfl_max = max_k |w_face(k)| * dt / h_eff(k),
  !  with h_eff(k) = 0.5*(dz(k)+dz(k+1)) for k=1..N-1.
  !
  !  nsubsteps is chosen so each substep has CFL_sub <= 1 approximately:
  !     nsubsteps = max(1, int(cfl_max) + 1)
  !--------------------------------------------------------------------
  pure subroutine compute_transport_substeps(w_face, dz, dt, cfl_max, nsubsteps)
      real(rk), intent(in)  :: w_face(0:,:)   ! (0:N, ntr)
      real(rk), intent(in)  :: dz(:)          ! (1:N)
      real(rk), intent(in)  :: dt
      real(rk), intent(out) :: cfl_max
      integer,  intent(out) :: nsubsteps

      real(rk), parameter :: cfl_target = 0.9_rk   ! target CFL per substep

      integer  :: N, ntr, k, ivar
      real(rk) :: h_eff, cfl_local, nsub_real

      N   = size(dz)
      ntr = size(w_face, 2)

      cfl_max = 0._rk
      nsubsteps = 1

      if (dt <= 0._rk .or. N <= 1 .or. ntr < 1) return
      if (size(w_face,1) < N+1) return

      ! Computing the Courant stability condition for each layer and tracer
      do k = 1, N-1
         h_eff = 0.5_rk * (dz(k) + dz(k+1))
         if (h_eff > 0._rk) then
            do ivar = 1, ntr
               cfl_local = abs(w_face(k, ivar)) * dt / h_eff
               if (cfl_local > cfl_max) cfl_max = cfl_local
            end do
         end if
      end do

      if (cfl_max <= 0._rk) then
        nsubsteps = 1
      else
        nsub_real = cfl_max / cfl_target   ! desired CFL_sub ~ cfl_target

        if (nsub_real <= 1._rk) then
            ! CFL already below target 
            nsubsteps = 1
        else
            nsubsteps = int(nsub_real) + 1
            if (nsubsteps < 1) nsubsteps = 1
        end if
      end if
   end subroutine compute_transport_substeps

  ! Compute the number of substeps needed for vertical mixing depending on the cnpar value
  pure subroutine compute_diffusion_substeps(Kz, dz, dt, cnpar, nsub, dt_sub)
    real(rk), intent(in)  :: Kz(0:), dz(:), dt, cnpar
    integer,  intent(out) :: nsub
    real(rk), intent(out) :: dt_sub

    integer  :: N, k
    real(rk) :: h_eff, Kloc, cflD, cflD_max

    N = size(dz)
    cflD_max = 0._rk

    do k = 1, N-1
        h_eff = 0.5_rk*(dz(k) + dz(k+1))
        Kloc  = max(Kz(k), 0._rk)

        if (h_eff > 0._rk .and. Kloc > 0._rk) then
            ! Diffusion CFL = (explicit fraction)*(dt)*(2K / h_eff^2)
            cflD = (1._rk - cnpar)*dt * 2._rk*Kloc / (h_eff*h_eff)
            if (cflD > cflD_max) cflD_max = cflD
        end if
    end do

    if (cflD_max <= 1._rk) then
        nsub = 1
        dt_sub = dt
    else
        ! We want cflD_sub = cflD_max / nsub ≤ 1, so nsub ≥ cflD_max
        nsub = int(cflD_max) + 1
        if (nsub < 1) nsub = 1
        dt_sub = dt / real(nsub, rk)
    end if
  end subroutine compute_diffusion_substeps


  !--------------------------------------------------------------------
  ! Compute number of bio substeps for explicit Euler integration
  ! based on a maximum change allowed per time-step.
  ! The idea is to keep integration using Forward Euler stable. 
  !
  !  - tendency_int(k,ivar): Tendencies for interor tracers [dC/dt]
  !  - tendency_sf(ivar):    tendencies for surface-only tracers [dc/dt]
  !  - tendency_bt(ivar):    tendencies for bottom-only tracers [dC/dt]
  !
  !  We estimate, for each variable:
  !      r = |dt_main * sms| / max(|C|, C_min)
  !  and choose nsub so that per substep r/nsub <= frac_max.
  !--------------------------------------------------------------------
  subroutine compute_bioint_substeps(BE, tendency_int, tendency_sf, tendency_bt, dt_main,          &
                                  frac_max, nsub, dt_sub)
      use bio_types,       only: BioEnv

      type(BioEnv), intent(in)  :: BE
      real(rk),     intent(in)  :: tendency_int(:,:), tendency_sf(:), tendency_bt(:)
      real(rk),     intent(in)  :: dt_main
      real(rk),     intent(in)  :: frac_max      ! e.g. 0.25_rk
      integer,      intent(out) :: nsub
      real(rk),     intent(out) :: dt_sub

      real(rk), parameter :: C_min   = 1.0e-30_rk   ! Just to avoid dividing by zero

      integer  :: nz, nint, nsfc, nbtm
      integer  :: k, ivar
      real(rk) :: C0, dC, r, r_max, nsub_real

      nz   = BE%grid%nz
      nint = BE%BS%n_interior
      nsfc = BE%BS%n_surface
      nbtm = BE%BS%n_bottom

      r_max = 0.0_rk

      ! Interuor variables
      if (nint>0) then
        do ivar = 1, nint
          do k = 1, nz
              C0   = BE%BS%interior_state(k,ivar)   ! Initial state
              dC   = dt_main * tendency_int(k,ivar) ! Change in the variable over the main time-step
              r    = abs(dC) / max(abs(C0), C_min)  ! How much the reported tendency is changing the state of the variable
              if (r > r_max) r_max = r              ! Maximum fractional change found for interior variables
          end do
        end do
      end if

      ! Surface-only variables
      if (nsfc>0) then
        do ivar = 1, nsfc
          C0   = BE%BS%surface_state(ivar)
          dC   = dt_main * tendency_sf(ivar)
          r    = abs(dC) / max(abs(C0), C_min)
          if (r > r_max) r_max = r
        end do
      end if

      ! Bottom-only tracers
      if (nbtm>0) then
        do ivar = 1, nbtm
          C0   = BE%BS%bottom_state(ivar)
          dC   = dt_main * tendency_bt(ivar)
          r    = abs(dC) / max(abs(C0), C_min)
          if (r > r_max) r_max = r
        end do
      end if

      ! Calculate nsub and dt_sub
      if (r_max <= frac_max .or. r_max <= 0._rk) then
        nsub = 1
      else
        nsub_real = r_max / frac_max
        nsub = ceiling(nsub_real)
        if (nsub < 1) nsub = 1
      end if

      dt_sub = dt_main / real(nsub, rk)

  end subroutine compute_bioint_substeps 


  !====================================================
  ! For Sediment processes
  !====================================================


  pure subroutine compute_bioirrigation_substeps(alpha, porewat_thickness, dz_wat_btm, dt_main, nsub)
      real(rk), intent(in)  :: alpha(:)              ! (nsed)
      real(rk), intent(in)  :: porewat_thickness(:)  ! (nsed) [m]
      real(rk), intent(in)  :: dz_wat_btm            ! [m]
      real(rk), intent(in)  :: dt_main               ! [s]
      integer,  intent(out) :: nsub

      real(rk), parameter :: r_target      = 0.9_rk
      real(rk), parameter :: lambda_target = 0.9_rk
      real(rk), parameter :: eps = 1.0e-30_rk

      real(rk) :: rmax, lambda, req
      real(rk) :: sum_exch

      rmax = maxval(max(alpha, 0._rk)) * dt_main

      sum_exch = sum(max(alpha, 0._rk) * max(porewat_thickness, 0._rk)) ! [m/s]
      lambda = dt_main * sum_exch / max(dz_wat_btm, eps)

      req = max(rmax / r_target, lambda / lambda_target)

      if (req <= 1._rk) then
        nsub = 1
      else
        nsub = ceiling(req)
      end if
  end subroutine compute_bioirrigation_substeps


  pure subroutine compute_deposition_substeps(BE, kwb, dz_w_bot, dt_main, nsub)
      type(BioEnv), intent(in) :: BE
      integer,      intent(in) :: kwb
      real(rk),     intent(in) :: dz_w_bot, dt_main
      integer,      intent(out):: nsub

      real(rk), parameter :: cfl_target = 0.9_rk
      real(rk), parameter :: eps = 1.0e-30_rk
      integer :: ivar, nint
      real(rk) :: vdep_max, dt_safe

      nint = BE%BS%n_interior
      vdep_max = 0._rk

      do ivar = 1, nint
        if (BE%tracer_info(ivar)%disable_transport) cycle
        if (.not. BE%tracer_info(ivar)%is_particulate) cycle

        ! vel_swi is negative downward
        vdep_max = max(vdep_max, max(0._rk, -BE%velocity(kwb, ivar)))
      end do

      if (vdep_max <= 0._rk) then
        nsub = 1
      else
        dt_safe = cfl_target * dz_w_bot / max(vdep_max, eps)
        nsub = max(1, ceiling(dt_main / dt_safe))
      end if
  end subroutine compute_deposition_substeps

  pure subroutine compute_diffusion_substeps_sed(diff, dz, phase_thickness, dt, cnpar, nsub, dt_sub)
      real(rk), intent(in)  :: diff(0:), dz(:), phase_thickness(:), dt, cnpar
      integer,  intent(out) :: nsub
      real(rk), intent(out) :: dt_sub

      integer :: N, i
      real(rk) :: a, c, cflD, cflD_max
      real(rk), parameter :: eps = 1.0e-30_rk

      N = size(dz)
      cflD_max = 0._rk

      if (N <= 1) then
        nsub = 1; dt_sub = dt; return
      end if
      if (cnpar >= 1._rk) then
        nsub = 1; dt_sub = dt; return
      end if

      do i = 2, N-1
        if (phase_thickness(i) <= eps) cycle
        a = 2._rk * max(diff(i-1),0._rk) / (dz(i)+dz(i-1)) / phase_thickness(i)
        c = 2._rk * max(diff(i  ),0._rk) / (dz(i)+dz(i+1)) / phase_thickness(i)
        cflD = (1._rk - cnpar) * dt * (a + c)
        if (cflD > cflD_max) cflD_max = cflD
      end do

      if (cflD_max <= 1._rk) then
        nsub = 1
      else
        nsub = int(cflD_max) + 1
        if (nsub < 1) nsub = 1
      end if

      dt_sub = dt / real(nsub, rk)
  end subroutine compute_diffusion_substeps_sed


  !--------------------------------------------------------------------
  ! Compute maximum vertical Courant number and numebr of substeps.
  ! This guarantees that each tracer never moves more than one layer each substep
  !
  !  - w_face(0:N) : vertical velocity at interfaces [m s-1]
  !  - dz(1:N)     : layer thickness [m]
  !  - dt          : time step [s]
  !
  !  cfl_max = max_k |w_face(k)| * dt / h_eff(k),
  !  with h_eff(k) = 0.5*(dz(k)+dz(k+1)) for k=1..N-1.
  !
  !  nsubsteps is chosen so each substep has CFL_sub <= 1 approximately:
  !     nsubsteps = max(1, int(cfl_max) + 1)
  !--------------------------------------------------------------------
  pure subroutine compute_transport_substeps_sed(w_face, dz, dt, cfl_max, nsubsteps)
      real(rk), intent(in)  :: w_face(0:)   ! (0:N)
      real(rk), intent(in)  :: dz(:)         ! (1:N)
      real(rk), intent(in)  :: dt
      real(rk), intent(out) :: cfl_max
      integer,  intent(out) :: nsubsteps

      real(rk), parameter :: cfl_target = 0.9_rk   ! target CFL per substep

      integer  :: N, k
      real(rk) :: h_eff, cfl_local, nsub_real

      N   = size(dz)
      cfl_max = 0._rk
      nsubsteps = 1

      if (dt <= 0._rk .or. N <= 1) return
      if (size(w_face) < N+1) return

      ! Computing the Courant stability condition for each layer and tracer
      do k = 1, N-1
         h_eff = 0.5_rk * (dz(k) + dz(k+1))
         if (h_eff > 0._rk) then            
            cfl_local = abs(w_face(k)) * dt / h_eff
            if (cfl_local > cfl_max) cfl_max = cfl_local
         end if
      end do

      if (cfl_max <= 0._rk) then
        nsubsteps = 1
      else
        nsub_real = cfl_max / cfl_target   ! desired CFL_sub ~ cfl_target

        if (nsub_real <= 1._rk) then
            ! CFL already below target 
            nsubsteps = 1
        else
            nsubsteps = max(1, ceiling(nsub_real))
            if (nsubsteps < 1) nsubsteps = 1
        end if
      end if
   end subroutine compute_transport_substeps_sed

end module numerical_stability
