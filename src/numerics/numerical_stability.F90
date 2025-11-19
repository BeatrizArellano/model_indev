module numerical_stability
  use precision_types, only: rk
  implicit none
  private
  public :: compute_phys_subcycles, compute_transport_substeps

  real(rk), parameter :: SAFETY      = 0.49_rk      ! inside parabolic limit with momentum diffusion explicit and Nz fixed over the subcycle
  real(rk), parameter :: DT_MIN      = 0.1_rk       ! [s] Stop below this time-step
  real(rk), parameter :: D_MIN       = 1.0e-12_rk   ! Minimum value for eddy coefficients so we don't divide by zero
  real(rk), parameter :: BOOST = 15._rk, FRAC_vismax = 0.7_rk

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
  end subroutine compute_phys_subcycles


  !--------------------------------------------------------------------
  ! Compute maximum vertical Courant number and numebr of substeps.
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
    real(rk), intent(in)  :: w_face(0:)
    real(rk), intent(in)  :: dz(:)
    real(rk), intent(in)  :: dt
    real(rk), intent(out) :: cfl_max
    integer,  intent(out) :: nsubsteps

    integer  :: N, k
    real(rk) :: h_eff, cfl_local

    N = size(dz)

    cfl_max = 0._rk

    ! Computing the Courant stability condition for each layer
    do k = 1, N-1
      ! Effective thickness around each interface
       h_eff = 0.5_rk * (dz(k) + dz(k+1))
       if (h_eff > 0._rk) then
          ! How many layers are crossed in one time-step
          cfl_local = abs(w_face(k)) * dt / h_eff
          if (cfl_local > cfl_max) cfl_max = cfl_local  ! Maximum Courant number over all layers
       end if
    end do

    if (cfl_max <= 0._rk) then
       nsubsteps = 1
    else
       nsubsteps = max(1, int(cfl_max) + 1)     ! Number of substeps needed to keep transport stable
       if (nsubsteps < 1) nsubsteps = 1
    end if

  end subroutine compute_transport_substeps

 

end module numerical_stability
