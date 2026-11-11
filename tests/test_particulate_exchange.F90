! Standalone numerical tests; compile with precision_types and sediment_exchange.
program test_particulate_exchange
    use precision_types,   only: rk
    use sediment_exchange, only: compute_natural_erosion, apply_particulate_exchange
    implicit none

    call test_erosion_law()
    call test_no_exchange()
    call test_deposition_only()
    call test_resuspension_only()
    call test_bidirectional()
    call test_zero_donors()
    call test_large_timestep()
    call test_substepping_consistency()
    call test_small_timestep()
    print '(A)', 'PASS: all sediment particulate exchange tests'
contains
    subroutine assert_close(actual, expected, label, rtol, atol)
        real(rk), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(rk), intent(in), optional :: rtol, atol
        real(rk) :: rel_tol, abs_tol

        rel_tol = 1.e-11_rk
        abs_tol = 1.e-13_rk
        if (present(rtol)) rel_tol = rtol
        if (present(atol)) abs_tol = atol
        if (abs(actual-expected) > abs_tol + rel_tol*max(abs(expected),abs(actual))) then
            print '(A,1X,A,2(1X,ES24.16))', 'FAIL:', label, actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_nonnegative(value, label)
        real(rk), intent(in) :: value
        character(*), intent(in) :: label
        if (value < 0._rk) then
            print '(A,1X,A,1X,ES24.16)', 'FAIL: negative', label, value
            error stop 1
        end if
    end subroutine assert_nonnegative

    subroutine test_erosion_law()
        real(rk) :: e
        call compute_natural_erosion(0.2_rk, 0.3_rk, 1.e-5_rk, e)
        call assert_close(e, 0._rk, 'stress below threshold')
        call compute_natural_erosion(0.3_rk, 0.3_rk, 1.e-5_rk, e)
        call assert_close(e, 0._rk, 'stress at threshold')
        call compute_natural_erosion(0.6_rk, 0.3_rk, 1.e-5_rk, e)
        call assert_close(e, 1.e-5_rk, 'stress twice threshold')
        call compute_natural_erosion(0.9_rk, 0.3_rk, 1.e-5_rk, e)
        call assert_close(e, 2.e-5_rk, 'stress three times threshold')
    end subroutine test_erosion_law

    subroutine test_no_exchange()
        real(rk) :: cw, cs, fd, fe
        cw = 2._rk
        cs = 10._rk
        call apply_particulate_exchange(cw, 2._rk, 0._rk, 3600._rk, 0.001_rk, &
            cs, 0.8_rk, 0._rk, 2650._rk, fd, fe)
        call assert_close(cw, 2._rk, 'no exchange water')
        call assert_close(cs, 10._rk, 'no exchange sediment')
        call assert_close(fd, 0._rk, 'no exchange deposition')
        call assert_close(fe, 0._rk, 'no exchange resuspension')
    end subroutine test_no_exchange

    subroutine test_deposition_only()
        real(rk) :: cw, cs, fd, fe, iw0, is0, iwnew, isnew, transferred
        real(rk), parameter :: dt=3600._rk, dw=2._rk, ds=0.001_rk, vel=-1.e-5_rk
        cw=2._rk
        cs=10._rk
        iw0=cw*dw
        is0=cs*ds
        call apply_particulate_exchange(cw,dw,vel,dt,ds,cs,0.8_rk,0._rk,2650._rk,fd,fe)
        iwnew=iw0*exp(vel*dt/dw)
        transferred=iw0-iwnew
        isnew=is0+transferred
        call assert_close(cw*dw,iwnew,'deposit-only water inventory')
        call assert_close(cs*ds,isnew,'deposit-only solid inventory')
        call assert_close(fd*dt,transferred,'deposit-only integrated flux')
        call assert_close(fe,0._rk,'deposit-only upward flux')
        call assert_close(cw*dw+cs*ds,iw0+is0,'deposit-only conserved')
    end subroutine test_deposition_only

    subroutine test_resuspension_only()
        real(rk) :: cw, cs, fd, fe, iw0, is0, isnew, transferred, ke
        real(rk), parameter :: dt=3600._rk, dw=2._rk, ds=0.001_rk, erosion=1.e-5_rk
        cw=2._rk
        cs=10._rk ! 50 tracer units / m3 solids when porosity is 0.8
        iw0=cw*dw
        is0=cs*ds
        ke=erosion/(2650._rk*0.2_rk*ds)
        call apply_particulate_exchange(cw,dw,0._rk,dt,ds,cs,0.8_rk,erosion,2650._rk,fd,fe)
        isnew=is0*exp(-ke*dt)
        transferred=is0-isnew
        call assert_close(cs*ds,isnew,'erosion-only sediment inventory')
        call assert_close(cw*dw,iw0+transferred,'erosion-only water inventory')
        call assert_close(fe*dt,transferred,'erosion-only integrated flux')
        call assert_close(fd,0._rk,'erosion-only deposition')
        call assert_close(cw*dw+cs*ds,iw0+is0,'erosion-only conserved')
    end subroutine test_resuspension_only

    subroutine test_bidirectional()
        real(rk) :: cw, cs, fd, fe, iw0, is0, kd, ke, k, iseq, expected_is
        real(rk), parameter :: dt=7200._rk, dw=2._rk, ds=0.001_rk, erosion=1.e-5_rk
        cw=0.5_rk
        cs=10._rk
        iw0=cw*dw
        is0=cs*ds
        kd=2.e-5_rk/dw
        ke=erosion/(2650._rk*0.2_rk*ds)
        k=kd+ke
        iseq=kd/k*(iw0+is0)
        expected_is=iseq+(is0-iseq)*exp(-k*dt)
        call apply_particulate_exchange(cw,dw,-2.e-5_rk,dt,ds,cs,0.8_rk,erosion,2650._rk,fd,fe)
        call assert_close(cs*ds,expected_is,'coupled analytic sediment')
        call assert_close(cw*dw+cs*ds,iw0+is0,'coupled conserved')
        call assert_close((cs*ds-is0)/dt,fd-fe,'coupled flux budget',atol=1.e-15_rk)
        call assert_nonnegative(cw,'coupled water')
        call assert_nonnegative(cs,'coupled sediment')
        if (fd <= 0._rk .or. fe <= 0._rk) error stop 'expected two positive gross fluxes'
    end subroutine test_bidirectional

    subroutine test_zero_donors()
        real(rk) :: cw, cs, fd, fe
        ! Newly deposited material may be resuspended within the same substep.
        cw=1._rk
        cs=0._rk
        call apply_particulate_exchange(cw,2._rk,-2.e-5_rk,3600._rk,0.001_rk,cs, &
            0.8_rk,1.e-5_rk,2650._rk,fd,fe)
        call assert_nonnegative(cw,'empty sediment initial water')
        call assert_nonnegative(cs,'empty sediment new concentration')
        if (fd <= 0._rk .or. fe <= 0._rk) error stop 'newly deposited material must recirculate'
        ! Material eroded earlier in the substep may be deposited again.
        cw=0._rk
        cs=10._rk
        call apply_particulate_exchange(cw,2._rk,-2.e-5_rk,3600._rk,0.001_rk,cs, &
            0.8_rk,1.e-5_rk,2650._rk,fd,fe)
        call assert_nonnegative(cw,'empty water initial water')
        call assert_nonnegative(cs,'empty water sediment')
        if (fd <= 0._rk .or. fe <= 0._rk) error stop 'newly eroded material must recirculate'
    end subroutine test_zero_donors

    subroutine test_large_timestep()
        real(rk) :: cw, cs, fd, fe, total, kd, ke, iseq
        cw=2._rk
        cs=10._rk
        total=cw*2._rk+cs*0.001_rk
        kd=2.e-5_rk/2._rk
        ke=1.e-5_rk/(2650._rk*0.2_rk*0.001_rk)
        iseq=kd/(kd+ke)*total
        call apply_particulate_exchange(cw,2._rk,-2.e-5_rk,1.e7_rk,0.001_rk,cs, &
            0.8_rk,1.e-5_rk,2650._rk,fd,fe)
        call assert_close(cw*2._rk+cs*0.001_rk,total,'long substep mass')
        call assert_close(cs*0.001_rk,iseq,'long substep equilibrium')
        call assert_nonnegative(cw,'long substep water')
        call assert_nonnegative(cs,'long substep sediment')
        if (fd <= 0._rk .or. fe <= 0._rk) error stop 'long substep gross fluxes missing'
    end subroutine test_large_timestep

    subroutine test_substepping_consistency()
        integer :: n
        real(rk) :: cw_big, cs_big, cw_small, cs_small, fd, fe, fd_big, fe_big
        real(rk) :: total_dep, total_resusp, small_dt
        real(rk), parameter :: dt=7200._rk, dw=2._rk, ds=0.001_rk
        cw_big=0.5_rk
        cs_big=10._rk
        cw_small=cw_big
        cs_small=cs_big
        call apply_particulate_exchange(cw_big,dw,-2.e-5_rk,dt,ds,cs_big, &
            0.8_rk,1.e-5_rk,2650._rk,fd_big,fe_big)
        total_dep=0._rk
        total_resusp=0._rk
        small_dt=dt/1000._rk
        do n=1,1000
            call apply_particulate_exchange(cw_small,dw,-2.e-5_rk,small_dt,ds,cs_small, &
                0.8_rk,1.e-5_rk,2650._rk,fd,fe)
            total_dep=total_dep+fd*small_dt
            total_resusp=total_resusp+fe*small_dt
        end do
        call assert_close(cw_small,cw_big,'substep invariance water')
        call assert_close(cs_small,cs_big,'substep invariance sediment')
        call assert_close(total_dep,fd_big*dt,'substep invariance gross deposition')
        call assert_close(total_resusp,fe_big*dt,'substep invariance gross resuspension')
    end subroutine test_substepping_consistency

    subroutine test_small_timestep()
        real(rk) :: cw, cs, fd, fe, kdep, ker, iw0, is0
        cw=2._rk
        cs=10._rk
        iw0=cw*2._rk
        is0=cs*0.001_rk
        kdep=1.e-5_rk/2._rk
        ker=1.e-5_rk/(2650._rk*0.2_rk*0.001_rk)
        call apply_particulate_exchange(cw,2._rk,-1.e-5_rk,1.e-9_rk,0.001_rk,cs, &
            0.8_rk,1.e-5_rk,2650._rk,fd,fe)
        call assert_close(fd,kdep*iw0,'tiny dt deposition rate',rtol=1.e-11_rk,atol=1.e-17_rk)
        call assert_close(fe,ker*is0,'tiny dt erosion rate',rtol=1.e-11_rk,atol=1.e-17_rk)
        call assert_close(cw*2._rk+cs*0.001_rk,iw0+is0,'tiny dt conserved')
    end subroutine test_small_timestep
end program test_particulate_exchange
