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
module vertical_mixing
  use precision_types, only: rk
  use tridiagonal, only: TridiagCoeff, solve_tridiag, reset_tridiag
  implicit none
  private

  public :: scalar_diffusion


contains

    ! ------ Vertical mixing parameterised as vertical diffusion -----------------
    !
    ! Solves vertical diffusion of a scalar variable (Var) on a layered water+sediment
    ! column using an implicit / semi-implicit scheme.
    !
    ! Assumptions:
    !   - Var(1:N) are layer-centre values, 1 = bottom, N = surface.
    !   - h(1:N)   are layer thicknesses [m].
    !   - Kz(0:N)  are interface diffusivities [m2 s-1], bottom..top.
    !       * Only Kz(1..N-1) (internal faces) are actually used here.
    !   - Boundary conditions:
    !       * Zero diffusive flux (homogeneous Neumann) at top and bottom.
    !       * No surface/bottom source terms or fluxes here – those are handled
    !         separately in reaction/flux steps.
    !
    ! cnpar:
    !   - Implicitness parameter:
    !       0.0  = explicit,
    !       0.5  = Crank–Nicolson,
    !       1.0  = fully implicit.
    !
    ! tricoef:
    !   - Reusable tridiagonal workspace (coefficients).
    !
    ! ierr: Status (0 = OK)
    subroutine scalar_diffusion(Var, N, dt, h, Kz, cnpar, tricoef, ierr)
        
        real(rk),           intent(inout) :: Var(1:N)
        integer,            intent(in)    :: N
        real(rk),           intent(in)    :: dt, cnpar
        real(rk),           intent(in)    :: h(1:N)            ! layer thicknesses
        real(rk),           intent(in)    :: Kz(0:N)           ! At layer interfaces: 0..N
        type(TridiagCoeff), intent(inout) :: tricoef
        integer,            intent(out)   :: ierr
  

        integer :: i
        real(rk) :: a, c
        
        ierr  = 0
        call reset_tridiag(tricoef)     

        !---------------------------
        ! Interior mixing: i = 2..N-1
        !---------------------------
        do i = 2, N-1
            c = 2.0_rk*dt*Kz(i)/(h(i)+h(i+1))/h(i)     ! couples to i+1
            a = 2.0_rk*dt*Kz(i-1)/(h(i)+h(i-1))/h(i)   ! couples to i-1

            tricoef%cu(i) = -cnpar * c
            tricoef%au(i) = -cnpar * a

            tricoef%bu(i) = 1.0_rk - (tricoef%au(i) + tricoef%cu(i))

            tricoef%du(i) = Var(i) + (1.0_rk-cnpar) * (a*Var(i-1) - (a+c)*Var(i) + c*Var(i+1))
        end do

        !---------------------------
        ! Bottom boundary: i = 1
        ! Neumann with zero diffusive flux
        !---------------------------
        c = 2.0_rk * dt * Kz(1) / (h(1) + h(2)) / h(1)

        tricoef%au(1) = 0.0_rk                      ! no layer below bottom
        tricoef%cu(1) = -cnpar * c                  ! couples to Var(2)
        tricoef%bu(1) = 1.0_rk - tricoef%cu(1)
        tricoef%du(1) = Var(1) + (1.0_rk - cnpar) * c * (Var(2) - Var(1))

        !---------------------------
        ! Top boundary: i = N
        ! Neumann with zero diffusive flux
        !---------------------------
        a = 2.0_rk * dt * Kz(N-1) / (h(N) + h(N-1)) / h(N)

        tricoef%cu(N) = 0.0_rk             ! no layer above surface
        tricoef%au(N) = -cnpar * a         ! couples to Var(N-1)
        tricoef%bu(N) = 1.0_rk - tricoef%au(N)

        tricoef%du(N) = Var(N) + (1.0_rk - cnpar) * a * (Var(N-1) - Var(N))     

        !---------------------------
        ! Solve tridiagonal system
        !---------------------------
        call solve_tridiag(1,N,tricoef,Var)

    end subroutine scalar_diffusion


end module vertical_mixing
