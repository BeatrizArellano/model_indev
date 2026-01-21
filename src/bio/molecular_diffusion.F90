module molecular_diffusion
  use precision_types, only: rk
  use bio_types,       only: TracerProperties, DIFF_NONE, DIFF_O2CO2_AB, DIFF_ION_LINEAR, &
                             DIFF_ARRHENIUS, DIFF_WILKE_CHANG, DIFF_STOKES_EINSTEIN 
  implicit none
  private

  public :: viscosity
  public :: molecular_diffusivity, run_diffusivity_test

  real(rk), parameter :: Patm_std = 1.013253_rk      ! [bar]
  real(rk), parameter :: R_gas    = 8.314462618_rk   ! [J mol^-1 K^-1]

contains

  !--------------------------------------------------------------------
  !> Dynamic viscosity of seawater [cP] using Boudreau (1997) Eq. 4.4
  !! Inputs:
  !!   t_C   - temperature [deg C]
  !!   S     - salinity [psu] (>=0)
  !!   P     - pressure [bar]
  !! Output:
  !!   mu_cP - dynamic viscosity [centipoise (i.e. 10-2 g cm-1 s-)]
  pure real(rk) function viscosity(t_C, S, P) result(mu_cP)
    real(rk), intent(in) :: t_C, S, P

    if (S < 0._rk) then
      mu_cP = 0._rk
      return
    end if

    mu_cP =  1.7910_rk - t_C*(6.144e-02_rk - t_C*(1.4510e-03_rk - t_C*1.6826e-05_rk)) &
           - 1.5290e-04_rk*P + 8.3885e-08_rk*P*P + 2.4727e-03_rk*S &
           + (6.0574e-06_rk*P - 2.6760e-09_rk*P*P)*t_C &
           + (t_C*(4.8429e-05_rk - t_C*(4.7172e-06_rk - t_C*7.5986e-08_rk)))*S

  end function viscosity


  !--------------------------------------------------------------------
  !> Molecular diffusivity in solution [m^2 s^-1] at in-situ conditions.
  !! Computes D(T,S,P) according to tr%diff_method.
  pure real(rk) function molecular_diffusivity(tr, t_C, S, P) result(D)
    type(TracerProperties), intent(in) :: tr
    real(rk), intent(in) :: t_C, S, P

    real(rk) :: TK, mu, mu0, mu_ref, mu0_poise, TrefK, D0

    D  = 0._rk
    TK = t_C + 273.15_rk                      ! Convert temperature from °C to °K
    mu = viscosity(t_C, S, P)                 ! [cP]
    mu0 = viscosity(t_C, 0._rk, Patm_std)     ! [cP] pure water @ 1 atm (needed in some methods)

    if (.not. tr%is_solute) return

    select case (tr%diff_method)

        case (DIFF_O2CO2_AB)
            ! Boudreau (1997) 4.58–4.59: D° = (A + B*(TK/mu0))*1e-5 [cm^2/s] then viscosity correction mu0/mu 
            D0 = (tr%A + tr%B * (TK / mu0)) * 1.0e-5_rk          ! [cm^2/s]
            D0 = D0 * 1.0e-4_rk                                  ! -> [m^2/s]
            D  = D0 * (mu0 / mu)                                 ! viscosity correction to (S,P)

        case (DIFF_ION_LINEAR)
            ! D° = (m0 + m1*t_C) * 1e-6 [cm^2/s]
            ! Coefficients m0, m1 are taken directly from Boudreau (1997) Tables 4.7 and 4.8. 
            D0 = (tr%m0 + tr%m1 * t_C) * 1.0e-6_rk               ! [cm^2/s]
            D0 = D0 * 1.0e-4_rk                                  ! -> [m^2/s]
            D  = D0 * (mu0 / mu)                                 ! viscosity correction to (S,P)

        case (DIFF_ARRHENIUS)
            ! D° = A0 * exp(-Ea/(R*TK)) [m^2/s]
            ! A0 in [1e-5 cm^2 s^-1], Ea in [kJ mol^-1]
            D0 = tr%A0 * exp(-(tr%Ea*1000._rk) / (R_gas * TK))   ! converts Ea from kJ/mol to J/mol [1e-5 cm^2/s]
            D0 = D0 * 1.0e-9_rk                                  ! -> [m^2/s]
            D = D0 * (mu0 / mu)                                  ! viscosity correction to (S,P)

        case (DIFF_WILKE_CHANG)
            ! Wilke–Chang (Wilke & Chang 1955; Hayduk & Laudie 1974 modification)
            ! 4.72E-09 * TK / (mu * Vb^0.6)  (Eq. 4.57 in Boudreau)
            ! Requires Vb [cm^3/mol], mu0 [in P]
            mu0_poise = mu0 * 1.0e-2_rk                          ! Convert from cP to Poise units
            D0 = 4.72e-09_rk * TK / (mu0_poise * tr%Vb**0.6_rk) 
            D0 = D0 * 1.0e-4_rk                                  ! -> [m^2/s]         
            D  = D0 * (mu0 / mu)                                 ! viscosity correction to (S,P)

        case (DIFF_STOKES_EINSTEIN)
            ! Uses a reference diffusivity and temperature to extrapolate to other temperatures
            ! D(T,S,P) = Dref * (TK/TrefK) * (mu(Tref,0,Patm) / mu(T,S,P))
            TrefK  = tr%Tref + 273.15_rk                         ! tr%Tref stored in [deg C]
            mu_ref = viscosity(tr%Tref, tr%Sref, Patm_std)
            D = tr%Dref * (TK / TrefK) * (mu_ref / mu)

        case default
            D = 0._rk
    end select

    if (D < 0._rk) D = 0._rk  
  end function molecular_diffusivity


  !========= TESTING subroutines ======================
  !> Create a TracerProperties instance for a named test tracer.
  !! This is a temporary check for validating the methods.
  pure function make_test_tracer(name) result(tr)
    character(len=*), intent(in) :: name
    type(TracerProperties) :: tr

    ! Defaults
    tr%fabm_index        = -1
    tr%is_solute         = .true.
    tr%is_particulate    = .false.
    tr%disable_transport = .false.
    tr%adsorption        = 0._rk

    tr%diff_method = 0
    tr%A  = 0._rk; tr%B  = 0._rk
    tr%m0 = 0._rk; tr%m1 = 0._rk
    tr%A0 = 0._rk; tr%Ea = 0._rk
    tr%Vb = 0._rk
    tr%Dref = 0._rk; tr%Tref = 0._rk; tr%Sref = 0._rk

    select case (trim(adjustl(name)))

    ! ---- O2 / CO2 (Boudreau Eqs 4.58-4.59)
    case ("O2")
      tr%diff_method = DIFF_O2CO2_AB
      tr%A = 0.2604_rk
      tr%B = 0.006383_rk

    case ("CO2")
      tr%diff_method = DIFF_O2CO2_AB
      tr%A = 0.1954_rk
      tr%B = 0.005089_rk

    ! ---- Ions (Boudreau Tables 4.7–4.8)  D°=(m0+m1*t)*1e-6 cm^2/s
    case ("NO3")
      tr%diff_method = DIFF_ION_LINEAR
      tr%m0 = 9.50_rk
      tr%m1 = 0.388_rk

    case ("SO4")
      tr%diff_method = DIFF_ION_LINEAR
      tr%m0 = 4.88_rk
      tr%m1 = 0.232_rk

    case ("NH4")
      tr%diff_method = DIFF_ION_LINEAR
      tr%m0 = 9.50_rk
      tr%m1 = 0.413_rk

    case ("HCO3")
      tr%diff_method = DIFF_ION_LINEAR
      tr%m0 = 5.06_rk
      tr%m1 = 0.275_rk

    ! ---- Arrhenius gases (Boudreau Table 4.4; Eq 4.60)
    ! A0 in [1e-5 cm^2/s], Ea in [kJ/mol]
    case ("H2")
      tr%diff_method = DIFF_ARRHENIUS
      tr%A0 = 3338._rk
      tr%Ea = 16.06_rk

    case ("CH4")
      tr%diff_method = DIFF_ARRHENIUS
      tr%A0 = 3047._rk
      tr%Ea = 18.36_rk

    case ("Ar")
      tr%diff_method = DIFF_ARRHENIUS
      tr%A0 = 7238._rk
      tr%Ea = 19.81_rk

    ! ---- Wilke–Chang (Boudreau Eq 4.57; Table 4.3)
    ! Vb in [cm^3/mol]
    case ("N2")
      tr%diff_method = DIFF_WILKE_CHANG
      tr%Vb = 34.7_rk

    case ("CO")
      tr%diff_method = DIFF_WILKE_CHANG
      tr%Vb = 34.5_rk

    case ("N2O")
      tr%diff_method = DIFF_WILKE_CHANG
      tr%Vb = 36.0_rk

    ! ---- Stokes–Einstein scaling (Boudreau Eq 4.107 form; marelac-style)
    ! Dref in [m^2/s], Tref [°C], Sref [psu]
    case ("H3PO4")
      tr%diff_method = DIFF_STOKES_EINSTEIN
      tr%Dref = 0.87e-9_rk
      tr%Tref = 25._rk
      tr%Sref = 0._rk

    case ("BOH3")
      tr%diff_method = DIFF_STOKES_EINSTEIN
      tr%Dref = 1.12e-9_rk
      tr%Tref = 25._rk
      tr%Sref = 29.2_rk

    case ("H4SiO4")
      tr%diff_method = DIFF_STOKES_EINSTEIN
      tr%Dref = 1.00e-9_rk
      tr%Tref = 25._rk
      tr%Sref = 36.1_rk

    case default
      ! Unknown tracer in this temporary test factory
      tr%is_solute    = .false.
      tr%diff_method  = 0
    end select
  end function make_test_tracer


  !> Simple test: prints diffusivities at a few canonical conditions.
    subroutine run_diffusivity_test()
        real(rk), parameter :: temps(*) = [0._rk, 10._rk, 25._rk]
        real(rk), parameter :: salts(*) = [0._rk, 35._rk]
        real(rk), parameter :: Pbar = Patm_std

        integer :: i, it, is
        type(TracerProperties) :: tr
        real(rk) :: D(size(temps), size(salts))
        logical :: ok_mono, ok_salt

        character(len=8), parameter :: TEST_TRACERS(*) = [ &
                                                            "O2     ", "CO2    ", &
                                                            "NO3    ", "SO4    ", "NH4    ", "HCO3   ", &
                                                            "H2     ", "CH4    ", "Ar     ", &
                                                            "N2     ", "CO     ", "N2O    ", &
                                                            "H3PO4  ", "BOH3   ", "H4SiO4 " ]


        write(*,'(a)') ""
        write(*,'(a)') "=============================================================="
        write(*,'(a)') "Molecular diffusivity smoke test "
        write(*,'(a)') "D returned in [m^2 s^-1]. Pressure for viscosity in [bar]."
        write(*,'(a,f10.6)') "Patm_std [bar] = ", Patm_std
        write(*,'(a)') "Temps [C] = 0,10,25 ; Salts = 0,35"
        write(*,'(a)') "=============================================================="

        do i = 1, size(TEST_TRACERS)
        tr = make_test_tracer(TEST_TRACERS(i))

        if (.not. tr%is_solute) then
            write(*,'(a,1x,a)') "SKIP (unknown tracer in factory?):", trim(TEST_TRACERS(i))
            cycle
        end if

        ! Compute matrix D(T,S)
        do is = 1, size(salts)
            do it = 1, size(temps)
            D(it,is) = molecular_diffusivity(tr, temps(it), salts(is), Pbar)
            end do
        end do

        ! Basic qualitative checks:
        ! (1) monotonic increase with T at each S
        ok_mono = .true.
        do is = 1, size(salts)
            if (D(2,is) < D(1,is)) ok_mono = .false.
            if (D(3,is) < D(2,is)) ok_mono = .false.
        end do

        ! (2) salinity effect at fixed T: D(S=35) should be <= D(S=0)
        ok_salt = .true.
        do it = 1, size(temps)
            if (D(it,2) > D(it,1)) ok_salt = .false.   ! salts(2)=35, salts(1)=0
        end do

        ! Print
        write(*,'(a)') ""
        write(*,'(a)') "Tracer: "//trim(TEST_TRACERS(i))
        write(*,'(a,i0)') "  diff_method = ", tr%diff_method
        if (.not. ok_mono) write(*,'(a)') "  WARNING: D not monotonic increasing with T (check units/viscosity)."
        if (.not. ok_salt) write(*,'(a)') "  WARNING: D at S=35 > S=0 for some T (check mu0/mu usage)."

        write(*,'(a)') "  Values (rows=T, cols=S):"
        write(*,'(a)') "            S=0               S=35"
        write(*,'(a,f5.1,2(1x,1pe14.6))') "  T= ", temps(1), D(1,1), D(1,2)
        write(*,'(a,f5.1,2(1x,1pe14.6))') "  T= ", temps(2), D(2,1), D(2,2)
        write(*,'(a,f5.1,2(1x,1pe14.6))') "  T= ", temps(3), D(3,1), D(3,2)

        end do

        write(*,'(a)') ""
        write(*,'(a)') "=============================================================="
        write(*,'(a)') "End of molecular diffusivity test."
        write(*,'(a)') "Stopping now"
        write(*,'(a)') "=============================================================="
        stop 0
    end subroutine run_diffusivity_test


end module molecular_diffusion
