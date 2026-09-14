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
! boundary conditions and optional implicit relaxation toward a target profile.
!
!=======================================================================================
module vertical_mixing
  use precision_types, only: rk
  use tridiagonal, only: TridiagCoeff, solve_tridiag, reset_tridiag
  implicit none
  private

  public :: scalar_diffusion
  public :: BC_DIRICHLET, BC_NEUMANN

  ! Types of boundary conditions for vertical diffusion
  integer, parameter :: BC_DIRICHLET = 1   ! prescribe the tracer value in the boundary cell
  integer, parameter :: BC_NEUMANN   = 2   ! prescribe flux into boundary cell [tracer m^-2 s^-1]

contains

    ! ------ Vertical mixing parameterised as vertical diffusion -----------------
    !
    ! Solves vertical diffusion of a scalar variable (Var) on a layered water column
    ! using an implicit / semi-implicit theta scheme.
    !
    ! Assumptions:
    !   - Var(1:N) are layer-centre values, 1 = bottom, N = surface.
    !   - h(1:N) are layer thicknesses [m].
    !   - diff(0:N) are interface diffusivities [m2 s-1], bottom..top.
    !     Only diff(1:N-1) (internal faces) are used by the diffusion operator.
    !
    ! cnpar:
    !   - Implicitness parameter for vertical diffusion only:
    !       0.0 = explicit,
    !       0.5 = Crank-Nicolson,
    !       1.0 = fully implicit.
    !
    ! Boundary conditions:
    !   - BC_NEUMANN: bc_*_value is flux INTO the boundary cell.
    !       Positive = source to the column, negative = sink from the column.
    !   - BC_DIRICHLET: bc_*_value is the prescribed boundary-cell concentration.
    !
    ! enforce_nonneg:
    !   If true, outward Neumann fluxes are linearised following Patankar (1980).
    !   The optional bc_*_flux_applied outputs report the realised Neumann flux
    !   after this linearisation. For non-Patankar Neumann boundaries they equal
    !   the requested flux; for Dirichlet boundaries they are returned as zero.
    !
    ! Relaxation:
    !   If relax_target and relax_timescale are both present, solve
    !
    !       dVar/dt = diffusion + (relax_target - Var) / relax_timescale
    !
    !   with the relaxation term treated fully implicitly (backward Euler),
    !   independently of cnpar. This contributes
    !
    !       bu += dt / relax_timescale
    !       du += (dt / relax_timescale) * relax_target
    !
    !   to each non-Dirichlet row. Hard Dirichlet boundary rows are not relaxed.
    !
    ! Defaults:
    !   - bc_top_type = BC_NEUMANN, bc_top_value = 0.0
    !   - bc_bot_type = BC_NEUMANN, bc_bot_value = 0.0
    subroutine scalar_diffusion(Var, N, dt, h, diff, cnpar, tricoef, &
                                ierr, enforce_nonneg,                &
                                bc_top_type, bc_top_value,           &
                                bc_bot_type, bc_bot_value,           &
                                relax_target, relax_timescale,       &
                                bc_top_flux_applied, bc_bot_flux_applied)

        real(rk),           intent(inout) :: Var(1:N)
        integer,            intent(in)    :: N
        real(rk),           intent(in)    :: dt, cnpar
        real(rk),           intent(in)    :: h(1:N)
        real(rk),           intent(in)    :: diff(0:N)
        type(TridiagCoeff), intent(inout) :: tricoef
        integer,            intent(out)   :: ierr

        ! Optionals
        logical,  intent(in),  optional :: enforce_nonneg
        integer,  intent(in),  optional :: bc_top_type, bc_bot_type
        real(rk), intent(in),  optional :: bc_top_value, bc_bot_value
        real(rk), intent(in),  optional :: relax_target(:)
        real(rk), intent(in),  optional :: relax_timescale
        real(rk), intent(out), optional :: bc_top_flux_applied, bc_bot_flux_applied

        integer :: i
        real(rk) :: a, c, relax_coeff
        real(rk), parameter :: tinyV = 1.0e-20_rk

        ! --- Local defaults/state ---
        logical  :: do_nonneg, do_relax
        logical  :: top_patankar, bot_patankar
        integer  :: top_type, bot_type
        real(rk) :: top_val, bot_val
        real(rk) :: top_denom, bot_denom

        do_nonneg = .false.
        do_relax  = .false.
        top_type  = BC_NEUMANN; top_val = 0.0_rk
        bot_type  = BC_NEUMANN; bot_val = 0.0_rk

        if (present(enforce_nonneg)) do_nonneg = enforce_nonneg
        if (present(bc_top_type))    top_type = bc_top_type
        if (present(bc_top_value))   top_val = bc_top_value
        if (present(bc_bot_type))    bot_type = bc_bot_type
        if (present(bc_bot_value))   bot_val = bc_bot_value

        if (present(bc_top_flux_applied)) bc_top_flux_applied = 0.0_rk
        if (present(bc_bot_flux_applied)) bc_bot_flux_applied = 0.0_rk

        ierr = 0

        ! Relaxation arguments must be supplied together.
        if (present(relax_target) .neqv. present(relax_timescale)) then
            ierr = 3
            return
        end if

        if (present(relax_target)) then
            if (size(relax_target) /= N) then
                ierr = 4
                return
            end if
            if (relax_timescale <= 0.0_rk) then
                ierr = 5
                return
            end if
            do_relax = .true.
            relax_coeff = dt / relax_timescale
        else
            relax_coeff = 0.0_rk
        end if

        top_patankar = do_nonneg .and. top_type == BC_NEUMANN .and. top_val < 0.0_rk
        bot_patankar = do_nonneg .and. bot_type == BC_NEUMANN .and. bot_val < 0.0_rk

        ! Store exactly the denominators used by the Patankar linearisation so
        ! the realised flux can be reconstructed consistently after the solve.
        top_denom = max(Var(N), tinyV)
        bot_denom = max(Var(1), tinyV)

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
            tricoef%du(i) = Var(i) + (1.0_rk-cnpar) * &
                            (a*Var(i-1) - (a+c)*Var(i) + c*Var(i+1))
        end do

        !---------------------------
        ! Bottom boundary: i = 1
        !---------------------------
        select case (bot_type)
        case (BC_NEUMANN)
            ! bot_val is flux INTO the bottom layer [Var m^-2 s^-1].
            c = 2.0_rk*dt*diff(1)/(h(1)+h(2))/h(1)
            tricoef%cu(1) = -cnpar * c

            if (bot_patankar) then
                tricoef%bu(1) = 1.0_rk - tricoef%cu(1) - dt*bot_val/(bot_denom*h(1))
                tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1))
            else
                tricoef%bu(1) = 1.0_rk - tricoef%cu(1)
                tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1)) + dt*bot_val/h(1)
            end if

        case (BC_DIRICHLET)
            tricoef%cu(1) = 0.0_rk
            tricoef%bu(1) = 1.0_rk
            tricoef%du(1) = bot_val

        case default
            ierr = 2
            return
        end select

        !---------------------------
        ! Top boundary: i = N
        !---------------------------
        select case (top_type)
        case (BC_NEUMANN)
            ! top_val is flux INTO the top layer [Var m^-2 s^-1].
            a = 2.0_rk*dt*diff(N-1)/(h(N)+h(N-1))/h(N)
            tricoef%au(N) = -cnpar * a

            if (top_patankar) then
                tricoef%bu(N) = 1.0_rk - tricoef%au(N) - dt*top_val/(top_denom*h(N))
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
            ierr = 1
            return
        end select

        !-------------------------------------------------------
        ! Optional relaxation to a prescribed water-column target.
        ! Fully implicit, following the GOTM diff_center formulation.
        ! Dirichlet rows remain hard prescribed values.
        !-------------------------------------------------------
        if (do_relax) then
            do i = 1, N
                if (i == 1 .and. bot_type == BC_DIRICHLET) cycle
                if (i == N .and. top_type == BC_DIRICHLET) cycle

                tricoef%bu(i) = tricoef%bu(i) + relax_coeff
                tricoef%du(i) = tricoef%du(i) + relax_coeff*relax_target(i)
            end do
        end if

        !---------------------------
        ! Solve tridiagonal system
        !---------------------------
        call solve_tridiag(1, N, tricoef, Var)

        !-----------------------------------------------------------------
        ! Report the Neumann flux actually represented by the solved matrix.
        ! For a Patankar sink the boundary term is proportional to the new
        ! concentration; otherwise the requested Neumann flux is unchanged.
        !-----------------------------------------------------------------
        if (present(bc_bot_flux_applied)) then
            if (bot_type == BC_NEUMANN) then
                if (bot_patankar) then
                    bc_bot_flux_applied = bot_val * Var(1) / bot_denom
                else
                    bc_bot_flux_applied = bot_val
                end if
            end if
        end if

        if (present(bc_top_flux_applied)) then
            if (top_type == BC_NEUMANN) then
                if (top_patankar) then
                    bc_top_flux_applied = top_val * Var(N) / top_denom
                else
                    bc_top_flux_applied = top_val
                end if
            end if
        end if

    end subroutine scalar_diffusion

end module vertical_mixing
