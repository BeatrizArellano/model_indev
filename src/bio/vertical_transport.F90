module vertical_transport
  use precision_types,      only: rk
  use grids,                only: VerticalGrid
  use numerical_stability,  only: compute_transport_substeps
  implicit none
  private

  public :: apply_vertical_transport

contains

  !--------------------------------------------------------------------
  ! Apply vertical transport due to residual movement (sinking or passive floating).
  !
  !  C(1:N)         : tracer at layer centres [units of C]
  !  grid%nz        : number of layers (1=bottom … N=surface)
  !  grid%dz(1:N)   : layer thicknesses [m]
  !  vel_center(1:N): Vertical velocity at layers' centres [m s-1]
  !                   (negative = sinking / downward; positive = upward)
  !  dt             : time step [s]
  !
  ! Assumptions:
  !  - Bottom-first indexing (1=bottom, N=surface).
  !  - No advective flux across bottom or surface for this process
  !    (w_face(0)=w_face(N)=0).
  !
  ! Uses the Conservative advection equation:
  !     dC/dt + d(vC)/dz = 0
  !
  ! The advection equation is integrated using a First-order upwind scheme.
  ! Concentrations are updated using Forward Euler
  ! 
  !
  ! On exit, C is updated in place.
  !--------------------------------------------------------------------
  subroutine apply_vertical_transport(C, grid, vel_center, dt)
    real(rk),           intent(inout) :: C(:)
    type(VerticalGrid), intent(in)    :: grid
    real(rk),           intent(in)    :: vel_center(:)
    real(rk),           intent(in)    :: dt

    integer  :: N, k, nsubsteps, istep
    real(rk) :: dt_sub, cfl_max
    real(rk), allocatable :: w_face(:)  ! 0..N
    real(rk), allocatable :: F(:)       ! Advective flux through the interfaces 0..N
    real(rk), allocatable :: dC(:)      ! Tendency per layer (dC/dt) 1..N

    ! --- Basic checks -------------------------------------------------
    N = grid%nz

    if (size(C) /= N) then
       stop 'apply_vertical_transport: size(C) /= grid%nz'
    end if

    if (size(vel_center) /= N) then
       stop 'apply_vertical_transport: size(vel_center) /= grid%nz'
    end if

    if (dt <= 0._rk) return

    allocate(w_face(0:N))
    allocate(F(0:N))
    allocate(dC(N))

    ! Compute velocities at layers' interfaces
    call velocity_at_interfaces(vel_center, grid%dz, w_face)

    ! Computes the stability condition (CFL) and the number of substeps
    ! this guarantees that each tracer never moves more than one layer each substep
    call compute_transport_substeps(w_face, grid%dz, dt, cfl_max, nsubsteps)
    dt_sub = dt / real(nsubsteps, rk)

    ! Inner loop
    do istep = 1, nsubsteps

       ! Compute fluxes at interfaces
       F(0) = 0._rk      ! no-flux at the bottom
       F(N) = 0._rk      ! no-flux at the surface

       ! From bottom to top
       do k = 1, N-1
          if (w_face(k) > 0._rk) then
             ! Upward flow: taking value from below (k)
             F(k) = w_face(k) * C(k)
          else if (w_face(k) < 0._rk) then
             ! Downward flow: take value from above (k+1)
             F(k) = w_face(k) * C(k+1)
          else
             F(k) = 0._rk
          end if
       end do

       ! Tendency from flux divergence (First-order upwind scheme)
       ! Discretisation of the 1D conservative advection equation
       do k = 1, N        
          ! Advection equation: dC/dt = -dF/dz  
          ! dC(k)/dt = - (1/dz(k))*(F(k)-F(k-1))
          dC(k) = - (F(k) - F(k-1)) / grid%dz(k)
       end do

       ! Update tracer concentrations using Forward Euler
       C(:) = C(:) + dt_sub * dC(:)

    end do

    deallocate(w_face, F, dC)

  end subroutine apply_vertical_transport

  !--------------------------------------------------------------------
  ! Compute velocities at interfaces from layer-centred vertical velocities.
  !
  !  vel_center(1:N) : velocity at centres [m s-1] (neg down, pos up)
  !  dz(1:N)         : layer thickness [m]
  !  w_face(0:N)     : velocity at interfaces [m s-1]
  !
  ! Bottom and surface are set to zero for now (no flux).
  !--------------------------------------------------------------------
  pure subroutine velocity_at_interfaces(vel_center, dz, w_face)
    real(rk), intent(in)  :: vel_center(:)
    real(rk), intent(in)  :: dz(:)
    real(rk), intent(out) :: w_face(0:)

    integer  :: N, k
    real(rk) :: alpha

    N = size(vel_center)

    ! No-flux boundaries for residual movement
    w_face(0) = 0._rk
    w_face(N) = 0._rk

    ! Thickness-weighted average for interior interfaces
    do k = 1, N-1
       alpha      = dz(k+1) / (dz(k) + dz(k+1))
       w_face(k)  = alpha * vel_center(k) + (1._rk - alpha) * vel_center(k+1)
    end do

  end subroutine velocity_at_interfaces

end module vertical_transport
