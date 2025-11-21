module vertical_mixing
  use precision_types, only: rk
  use physics_types,   only: Dirichlet, Neumann
  use tridiagonal, only: TridiagCoeff, solve_tridiag, reset_tridiag
  implicit none
  private

  public :: scalar_diffusion


contains

    ! ------ Vertical mixing parameterised as vertical diffusion -----------------
    ! Solves vertical diffusion of a scalar (Y) on a layered water column.
    ! Uses an implicit / semi-implicit approach for stability.
    ! Accepts face-Kz and variable layer thickness; supports Neumann/Dirichlet BCs.
    ! Optional: positivity guard at Neumann BCs, linear/constant sources, relaxation.
    ! Defaults mimick S2P3: zero-flux Neumann at top and bottom, no sources/relaxation.
    ! Y(1:N)           : variable at layer centers to be mixed [units of Y]
    ! N                : number of layers (>=3); centers indexed 1..N
    ! dt               : time step [s]
    ! h(1:N)           : layer thicknesses [m]
    ! Kz(0:N)          : vertical diffusivity at layer interfaces (bottom..top) [m^2 s^-1]
    ! cnpar            : implicitness (0=explicit, 0.5=Crank–Nicolson, 1=implicit)
    ! tricoef          : reusable tridiagonal workspace (coefficients + scratch)
    ! ierr             : status flag (0=ok; >0=BC error code; <0=input/other)
    !
    ! bc_top_type      : boundary type at surface (Dirichlet or Neumann)
    ! bc_top_value     : if Neumann: flux into top cell [Y·m s^-1]; if Dirichlet: prescribed Y at surface [Y]
    ! bc_bot_type      : boundary type at bottom (Dirichlet or Neumann)
    ! bc_bot_value     : if Neumann: flux into bottom cell [Y·m s^-1]; if Dirichlet: prescribed Y at bottom [Y]
    !
    ! enforce_nonneg   : if true, apply Patankar approach at Neumann BCs to prevent Y<0 (for positive tracers)

    ! Yobs(1:N)(opt)   : relaxation target profile [Y]
    ! L(1:N)   (opt)   : linear source coefficient (implicit) [s^-1]
    ! Q(1:N)   (opt)   : constant source term (explicit) [Y s^-1]
    ! Taur(1:N)(opt)   : relaxation timescale toward Yobs [s]     


    ! bc_*_value meaning:
    !   Neumann: flux into the boundary cell [Y·m s^-1].
    !   Dirichlet: prescribed value of Y at the boundary cell center.
    subroutine scalar_diffusion(Var, N, dt, h, Kz, cnpar,  &
                                tricoef, ierr, enforce_nonneg, &
                                 bc_top_type, bc_top_value,     &
                                 bc_bot_type, bc_bot_value,     &
                                 Yobs, Taur, LinTerm, Q)
        
        real(rk),           intent(inout) :: Var(1:N)
        integer,            intent(in)    :: N
        real(rk),           intent(in)    :: dt, cnpar
        real(rk),           intent(in)    :: h(1:N)            ! layer thicknesses
        real(rk),           intent(in)    :: Kz(0:N)           ! faces: 0..N
        type(TridiagCoeff), intent(inout) :: tricoef
        integer,            intent(out)   :: ierr
        ! Optionals
        logical,         intent(in), optional :: enforce_nonneg
        integer,         intent(in), optional :: bc_top_type, bc_bot_type
        real(rk),        intent(in), optional :: bc_top_value, bc_bot_value
        real(rk),        intent(in), optional :: Yobs(1:N), Taur(1:N)
        real(rk),        intent(in), optional :: LinTerm(1:N), Q(1:N)
  

        integer :: i
        real(rk) :: a, c, l
        real(rk), parameter :: tinyY = 1.0e-20_rk

        ! --- Local defaults ---
        logical  :: do_nonneg
        integer  :: top_type, bot_type
        real(rk) :: top_val,  bot_val        
        logical  :: hasL, hasQ, has_relax

        do_nonneg = .false.
        top_type  = Neumann;  top_val = 0.0_rk
        bot_type  = Neumann;  bot_val = 0.0_rk
        ! Read optionals (if provided)
        if (present(enforce_nonneg)) do_nonneg = enforce_nonneg
        if (present(bc_top_type))    top_type  = bc_top_type
        if (present(bc_top_value))   top_val   = bc_top_value
        if (present(bc_bot_type))    bot_type  = bc_bot_type
        if (present(bc_bot_value))   bot_val   = bc_bot_value

        hasL      = present(LinTerm)
        hasQ      = present(Q)
        has_relax = present(Taur) .and. present(Yobs)

        ierr  = 0



        call reset_tridiag(tricoef)      
        

        !---------------------------
        ! Interior: i = 2..N-1
        !---------------------------
        do i = 2, N-1
            c = 2.0_rk*dt*Kz(i)/(h(i)+h(i+1))/h(i)     ! couples to i+1
            a = 2.0_rk*dt*Kz(i-1)/(h(i)+h(i-1))/h(i)   ! couples to i-1

            tricoef%cu(i) = -cnpar * c
            tricoef%au(i) = -cnpar * a

            l = 0.0_rk; if (hasL) l = dt*LinTerm(i)
            tricoef%bu(i) = 1.0_rk - (tricoef%au(i) + tricoef%cu(i)) - l

            tricoef%du(i) = Var(i) + (1.0_rk-cnpar) * ( a*Var(i-1) - (a+c)*Var(i) + c*Var(i+1) )
            if (hasQ) tricoef%du(i) = tricoef%du(i) + dt*Q(i)
        end do

        !---------------------------
        ! Top boundary: i = N
        !---------------------------
        select case (top_type)
            case (Neumann)
                a = 2.0_rk*dt*Kz(N-1)/(h(N)+h(N-1))/h(N)
                tricoef%au(N) = -cnpar * a
                ! Patankar tweak if flux leaves the column and non-negativity enforced
                if (do_nonneg .and. top_val < 0.0_rk) then
                    tricoef%bu(N) = 1.0_rk - tricoef%au(N) - dt*top_val / (max(Var(N),tinyY)*h(N))
                    tricoef%du(N) = Var(N) + (1.0_rk-cnpar)*a*(Var(N-1)-Var(N))
                else
                    tricoef%bu(N) = 1.0_rk - tricoef%au(N)
                    tricoef%du(N) = Var(N) + (1.0_rk-cnpar)*a*(Var(N-1)-Var(N)) + dt*top_val/h(N)
                end if
                if (hasL) then
                    tricoef%bu(N) = tricoef%bu(N) - dt*LinTerm(N)    ! implicit linear source
                end if
                if (hasQ) then
                    tricoef%du(N) = tricoef%du(N) + dt*Q(N)    ! explicit constant source
                end if

            case (Dirichlet)
                tricoef%au(N) = 0.0_rk
                tricoef%bu(N) = 1.0_rk
                tricoef%du(N) = top_val

            case default
                ierr = 1; return
        end select

        !---------------------------
        ! Bottom boundary: i = 1
        !---------------------------
        select case (bot_type)
            case (Neumann)
                c = 2.0_rk*dt*Kz(1)/(h(1)+h(2))/h(1)
                tricoef%cu(1) = -cnpar * c
                if (do_nonneg .and. bot_val < 0.0_rk) then
                    tricoef%bu(1) = 1.0_rk - tricoef%cu(1) - dt*bot_val / (max(Var(1),tinyY)*h(1))
                    tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1))
                else
                    tricoef%bu(1) = 1.0_rk - tricoef%cu(1)
                    tricoef%du(1) = Var(1) + (1.0_rk-cnpar)*c*(Var(2)-Var(1)) + dt*bot_val/h(1)
                end if
                if (hasL) then
                    tricoef%bu(1) = tricoef%bu(1) - dt*LinTerm(1)
                end if
                if (hasQ) then
                    tricoef%du(1) = tricoef%du(1) + dt*Q(1)
                end if

            case (Dirichlet)
                tricoef%cu(1) = 0.0_rk
                tricoef%bu(1) = 1.0_rk
                tricoef%du(1) = bot_val

            case default
            ierr = 2; return
        end select

        !---------------------------
        ! Optional relaxation
        !---------------------------
        if (has_relax) then
            do i = 1, N
                if (Taur(i) < 1.0e10_rk) then
                    tricoef%bu(i) = tricoef%bu(i) + dt / Taur(i)
                    tricoef%du(i) = tricoef%du(i) + dt / Taur(i) * Yobs(i)
                end if
            end do
        end if

        !---------------------------
        ! Solve tridiagonal system
        !---------------------------
        call solve_tridiag(1,N,tricoef,Var)

    end subroutine scalar_diffusion


end module vertical_mixing
