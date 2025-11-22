module vertical_transport
  use precision_types,      only: rk
  use grids,                only: VerticalGrid
  implicit none
  private

  public :: apply_vertical_transport, velocity_at_interfaces

contains

   !--------------------------------------------------------------------
   ! Apply vertical transport (residual movement) for a SINGLE tracer.
   !
   !  C(1:N)        : tracer at layer centres [units of C], updated in place
   !  grid%nz       : number of layers (1=bottom … N=surface)
   !  grid%dz(1:N)  : layer thicknesses [m]
   !  w_face(0:N)   : vertical velocity at interfaces [m s-1]
   !                  (neg = downward, pos = upward). 
   !  dt            : timestep
   !
   ! Assumptions:
   !  - Bottom-first indexing (1=bottom, N=surface).
   !  - No advective flux across bottom or surface for this process
   !    (w_face(0) and w_face(N) are assumed to be 0).
   !
   ! Uses the Conservative advection equation:
   !     dC/dt + d(vC)/dz = 0
   !
   ! The advection equation is integrated using a First-order upwind scheme.
   ! Concentrations are updated using Forward Euler
   !
   ! On exit, C is updated in place.
   !--------------------------------------------------------------------
   subroutine apply_vertical_transport(C, grid, w_face, dt)
      real(rk),           intent(inout) :: C(:)
      type(VerticalGrid), intent(in)    :: grid
      real(rk),           intent(in)    :: w_face(0:)
      real(rk),           intent(in)    :: dt

      integer  :: N, k
      real(rk), allocatable :: F(:)   ! fluxes at interfaces 0..N
      real(rk), allocatable :: dC(:)  ! tendencies 1..N

      N = grid%nz

      ! Basic checks
      if (size(C) /= N) then
         stop 'apply_vertical_transport: size(C) /= grid%nz'
      end if
      if (size(w_face) /= N+1) then
         stop 'apply_vertical_transport: size(w_face) /= grid%nz+1'
      end if
      if (dt <= 0._rk) return

      allocate(F(0:N))
      allocate(dC(N))

      
      F(0) = 0._rk    ! No flux at bottom
      F(N) = 0._rk    ! No flux at surface

      ! Compute fluxes at interior interfaces (first-order upwind)
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
         dC(k) = - (F(k) - F(k-1)) / grid%dz(k)
      end do

      ! Update tracer concentrations using Forward Euler
      C(:) = C(:) + dt * dC(:)

      deallocate(F, dC)

   end subroutine apply_vertical_transport
   


   !--------------------------------------------------------------------
   ! Compute velocities at interfaces from layer-centred vertical
   ! velocities, for ALL tracers at once.
   !
   !  vel_center(1:N,1:ntr) : velocity at centres [m s-1]
   !                          (neg = downward, pos = upward)
   !  grid%dz(1:N)          : layer thickness [m]
   !  w_face(0:N,1:ntr)     : velocity at interfaces [m s-1]
   !
   ! Bottom and surface are set to zero for now (no-flux residual movement).
   !--------------------------------------------------------------------
   pure subroutine velocity_at_interfaces(vel_center, grid, w_face)
      real(rk),           intent(in)  :: vel_center(:,:)   ! (1:N, 1:ntr)
      type(VerticalGrid), intent(in)  :: grid
      real(rk),           intent(out) :: w_face(0:,:)      ! (0:N, 1:ntr)

      integer  :: N, ntr, k, ivar
      real(rk) :: alpha

      N   = grid%nz
      ntr = size(vel_center, 2)

      ! No-flux boundaries for residual movement
      do ivar = 1, ntr
         w_face(0,ivar) = 0._rk
         w_face(N,ivar) = 0._rk
      end do

      ! Velocities at interfaces weighed by the layer thickness
      do ivar = 1, ntr
         do k = 1, N-1
            alpha        = grid%dz(k+1) / (grid%dz(k) + grid%dz(k+1))
            w_face(k,ivar) = alpha*vel_center(k,ivar) + (1._rk - alpha)*vel_center(k+1,ivar)
         end do
      end do

   end subroutine velocity_at_interfaces

end module vertical_transport
