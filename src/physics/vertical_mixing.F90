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
! scalar_diffusion is modified in this version to handle Neumann and Dirichlet 
!   boundary conditions. 
!
!=======================================================================================
module vertical_mixing
  use precision_types, only: rk
  use tridiagonal, only: TridiagCoeff, solve_tridiag, reset_tridiag
  implicit none
  private

  public :: scalar_diffusion
  public :: BC_DIRICHLET, BC_NEUMANN

  ! Types of boundary conditiosn for vertical diffusion
  integer, parameter :: BC_DIRICHLET = 1   ! BC_DIRICHLET: prescribe the tracer value at the boundary (fixed concentration)
  integer, parameter :: BC_NEUMANN = 2     ! BC_NEUMANN  : prescribe the tracer flux at the boundary (flux into cell is positive) [tracer m^-2 s^-1]

contains

    ! ------ Vertical mixing parameterised as vertical diffusion -----------------
    !
    ! Solves vertical diffusion of a scalar variable (Var) on a layered water+sediment
    ! column using an implicit / semi-implicit scheme.
    !
    ! Assumptions:
    !   - Var(1:N) are layer-centre values, 1 = bottom, N = surface.
    !   - h(1:N)   are layer thicknesses [m].
    !   - diff(0:N)  are interface diffusivities [m2 s-1], bottom..top.
    !       * Only diff(1..N-1) (internal faces) are actually used here.
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
    !
    ! Boundary conditions
    !
    ! bc_top_type, bc_bot_type:
    !   - BC_NEUMANN  : prescribe a diffusive flux at the boundary
    !   - BC_DIRICHLET: prescribe the tracer value in the boundary layer
    !
    ! bc_top_value, bc_bot_value meaning:
    !
    !   * For BC_NEUMANN:
    !       The value represents a flux INTO the boundary cell
    !       with units [Var m^-2 s^-1].
    !       - Positive value  -> source to the column
    !       - Negative value  -> sink from the column
    !
    !       The flux contribution is applied as:
    !           dVar/dt += flux / h(cell)
    !       i.e. the flux is converted to a volumetric tendency.
    !
    !   * For BC_DIRICHLET:
    !       The value represents the tracer concentration prescribed
    !       in the boundary cell (cell-centre value). This modifies
    !       the boundary layer value directly
    ! enforce_nonneg:
    !   If true, enforce non-negativity for tracers by applying
    !   the Patankar (1980) linearisation for outward Neumann fluxes.
    !
    !   For BC_NEUMANN with flux < 0 (flux leaving the column),
    !   the flux is treated as a linear sink proportional to the
    !   local concentration and added implicitly to the diagonal.
    !
    ! Defaults:
    !   - bc_top_type = BC_NEUMANN, bc_top_value = 0.0   (zero diffusive flux)
    !   - bc_bot_type = BC_NEUMANN, bc_bot_value = 0.0   (zero diffusive flux)
    subroutine scalar_diffusion(Var, N, dt, h, diff, cnpar, tricoef,  &
                                ierr, enforce_nonneg,                 &
                                bc_top_type, bc_top_value,            &
                                bc_bot_type, bc_bot_value)
        
        real(rk),           intent(inout) :: Var(1:N)
        integer,            intent(in)    :: N
        real(rk),           intent(in)    :: dt, cnpar
        real(rk),           intent(in)    :: h(1:N)              ! layer thicknesses
        real(rk),           intent(in)    :: diff(0:N)           ! At layer interfaces: 0..N
        type(TridiagCoeff), intent(inout) :: tricoef
        integer,            intent(out)   :: ierr

        ! Optionals
        logical,         intent(in), optional :: enforce_nonneg
        integer,         intent(in), optional :: bc_top_type, bc_bot_type
        real(rk),        intent(in), optional :: bc_top_value, bc_bot_value
  

        integer :: i
        real(rk) :: a, c
        real(rk), parameter :: tinyV = 1.0e-20_rk

        ! --- Local defaults ---
        logical  :: do_nonneg
        integer  :: top_type, bot_type
        real(rk) :: top_val,  bot_val  

        do_nonneg = .false.
        top_type  = BC_NEUMANN;  top_val = 0.0_rk
        bot_type  = BC_NEUMANN;  bot_val = 0.0_rk
        ! Read optionals (if provided)
        if (present(enforce_nonneg)) do_nonneg = enforce_nonneg
        if (present(bc_top_type))    top_type  = bc_top_type
        if (present(bc_top_value))   top_val   = bc_top_value
        if (present(bc_bot_type))    bot_type  = bc_bot_type
        if (present(bc_bot_value))   bot_val   = bc_bot_value
        
        ierr  = 0
        call reset_tridiag(tricoef)     

        !---------------------------
        ! Interior mixing: i = 2..N-1
        !---------------------------
        do i = 2, N-1
            c = 2.0_rk*dt*diff(i)/(h(i)+h(i+1))/h(i)     ! couples to i+1
            a = 2.0_rk*dt*diff(i-1)/(h(i)+h(i-1))/h(i)   ! couples to i-1

            tricoef%cu(i) = -cnpar * c
            tricoef%au(i) = -cnpar * a

            tricoef%bu(i) = 1.0_rk - (tricoef%au(i) + tricoef%cu(i))

            tricoef%du(i) = Var(i) + (1.0_rk-cnpar) * (a*Var(i-1) - (a+c)*Var(i) + c*Var(i+1))
        end do

        !---------------------------
        ! Bottom boundary: i = 1
        !---------------------------
        select case (bot_type)
            case (BC_NEUMANN)
                !   bot_val is the flux INTO the bottom layer [Var m^-2 s^-1]
                !   (positive = into the layer, negative = out of the domain)
                c = 2.0_rk* dt * diff(1)/(h(1)+h(2))/h(1)
                tricoef%cu(1) = -cnpar * c
                if (do_nonneg .and. bot_val < 0.0_rk) then
                    ! Patankar (1980) approach if flux leaves the column and non-negativity enforced
                    tricoef%bu(1) = 1.0_rk - tricoef%cu(1) - dt*bot_val / (max(Var(1),tinyV)*h(1))
                    tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1))
                else
                    tricoef%bu(1) = 1.0_rk - tricoef%cu(1)
                    tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1)) + dt*bot_val/h(1)
                end if

            case (BC_DIRICHLET)
                !   The boundary value is prescribed directly.
                !   This sets Var(1) / Var(N) to the specified value
                tricoef%cu(1) = 0.0_rk
                tricoef%bu(1) = 1.0_rk
                tricoef%du(1) = bot_val

            case default
            ierr = 2; return
        end select

        !---------------------------
        ! Top boundary: i = N
        !---------------------------
        select case (top_type)
            case (BC_NEUMANN)
                !   top_val is the flux INTO the top layer [Var m^-2 s^-1]
                !   (positive = into the layer, negative = out of the domain)
                a = 2.0_rk * dt * diff(N-1)/(h(N)+h(N-1))/h(N)
                tricoef%au(N) = -cnpar * a
                ! Patankar approach if flux leaves the column and non-negativity enforced
                if (do_nonneg .and. top_val < 0.0_rk) then
                    tricoef%bu(N) = 1.0_rk - tricoef%au(N) - dt*top_val / (max(Var(N),tinyV)*h(N))
                    tricoef%du(N) = Var(N) + (1.0_rk-cnpar)*a*(Var(N-1)-Var(N))
                else
                    tricoef%bu(N) = 1.0_rk - tricoef%au(N)
                    tricoef%du(N) = Var(N) + (1.0_rk-cnpar)*a*(Var(N-1)-Var(N)) + dt*top_val/h(N)
                end if

            case (BC_DIRICHLET)
                tricoef%au(N) = 0.0_rk
                tricoef%bu(N) = 1.0_rk
                tricoef%du(N) = top_val

            case default
                ierr = 1; return
        end select    
        
        !---------------------------
        ! Solve tridiagonal system
        !---------------------------
        call solve_tridiag(1,N,tricoef,Var)

    end subroutine scalar_diffusion


end module vertical_mixing
