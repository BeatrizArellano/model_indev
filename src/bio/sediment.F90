module sediment    
    use bio_types,         only: SedimentEnv, BioEnv, TracerProperties, DIFF_NONE, DIFF_O2CO2_AB, &
                                 DIFF_ION_LINEAR, DIFF_ARRHENIUS, DIFF_WILKE_CHANG, DIFF_STOKES_EINSTEIN
    use bio_params,        only: SedParams, read_sed_parameters
    use grids,             only: VerticalGrid
    use precision_types,   only: rk
    use read_config_yaml,  only: ConfigParams
    use output_static,     only: StaticProfile, init_static_profile, append_static_profile
    use str_utils,         only: to_lower
    use tridiagonal,       only: init_tridiag, clear_tridiag

    implicit none
    private

    public :: init_sediment, output_sed_profiles, clear_sediment_env, write_tracer_properties
    public :: phase_to_bulk, bulk_to_phase, phase_to_bulk_all, bulk_to_phase_all

contains

    ! Initialise sediment environment for one column.
    !! Reads user parameters from YAML, converts to SI units, allocates arrays,
    !! and precomputes sediment property profiles (porosity, tortuosity, Db, irrigation).
    subroutine init_sediment(cfg_params, grid, SE)
        type(ConfigParams),  intent(in) :: cfg_params
        type(VerticalGrid),  intent(in), target :: grid        ! Sediment grid
        type(SedimentEnv),   intent(inout) :: SE

        ! Associate grid pointer
        SE%grid => grid
        SE%nz = grid%nz

        ! Read sediment parameters
        call read_sed_parameters(cfg_params, SE%params_user)
        ! Convert parameter units to SI units cm->m and yr->s
        call convert_units_to_SI(SE%params_user, SE%params_SI)
        ! Allocate working arrays 
        call allocate_sed_arrays(SE)

        SE%use_bioturbation       = (SE%params_SI%biot_mode /= 'off')
        SE%use_bioirrigation      = (SE%params_SI%irr_mode  /= 'off')
        SE%output_bioturb_dynamic = (SE%params_SI%biot_mode == 'dynamic')
        SE%output_bioirr_dynamic  = (SE%params_SI%irr_mode  == 'dynamic')

        ! Compute profiles for sediment properties
        call compute_porosity_profile(grid, SE%params_SI, SE%poro, SE%poro_w, SE%theta2, SE%porewat_thickness, SE%solid_thickness)
        if (SE%use_bioturbation) then
            call compute_bioturbation_profile(grid, SE%params_SI%biot_db_sfc, SE%params_SI%biot_mld, SE%params_SI%biot_ez, SE%bioturb)
        else 
            SE%bioturb = 0.0_rk
        end if
        if (SE%use_bioirrigation) then
            call compute_bioirrigation_alpha(grid, SE%params_SI%irr_sfc, SE%params_SI%irr_ez, SE%bioirr, SE%bioirr_w)
        else 
            SE%bioirr   = 0.0_rk
            SE%bioirr_w = 0.0_rk
        end if
        ! ---- Burial velocities (interfaces). 
        call compute_burial_velocities(SE%grid, SE%params_SI, SE%poro_w, SE%params_SI%sed_rate, SE%vel_solids, SE%vel_solutes)

        ! Initialise tridiagonal matrix
        call init_tridiag(SE%sed_trid, SE%nz)

        SE%is_init = .true.
        ! call write_sediment_profiles('Sediment_profiles.dat',SE%grid, SE)
    end subroutine init_sediment

    ! Convert from phase-specific concentrations to bulk-sediment concentrations
    ! For solids (bulk=(1-phi)*C) and for solutes (bulk=(phi)*C)
    pure subroutine phase_to_bulk(C_phase, poro, is_solute, C_bulk)
        real(rk), intent(in)  :: C_phase(:)      ! phase-specific concentration (nsed)
        real(rk), intent(in)  :: poro(:)         ! porosity at centres (nsed)
        logical, intent(in)   :: is_solute
        real(rk), intent(out) :: C_bulk(:)       ! bulk-sediment concentration (nsed)

        integer :: n, k
        real(rk), parameter :: tiny = 1.0e-12_rk
        real(rk) :: phi, one_m_phi

        n = size(C_phase)

        if (is_solute) then
            do k = 1, n
                phi = max(poro(k), tiny)
                C_bulk(k) = phi * C_phase(k)
            end do
        else
            do k = 1, n
                one_m_phi = max(1.0_rk - poro(k), tiny)
                C_bulk(k) = one_m_phi * C_phase(k)
            end do
        end if
    end subroutine phase_to_bulk

    ! Convert from bulk-sediment concentrations to phase-specific concentrations
    ! For solids (bulk=(1-phi)*C) and for solutes (bulk=(phi)*C)
    pure subroutine bulk_to_phase(C_bulk, poro, is_solute, C_phase)
        real(rk), intent(in)  :: C_bulk(:)       ! bulk-sediment concentration (nsed)
        real(rk), intent(in)  :: poro(:)         ! porosity at centres (nsed)
        logical, intent(in)   :: is_solute
        real(rk), intent(out) :: C_phase(:)      ! phase-specific concentration (nsed)

        integer :: n, k
        real(rk), parameter :: tiny = 1.0e-12_rk
        real(rk) :: phi, one_m_phi

        n = size(C_bulk)

        if (is_solute) then
            do k = 1, n
                phi = max(poro(k), tiny)
                C_phase(k) = C_bulk(k) / phi
            end do
        else
            do k = 1, n
                one_m_phi = max(1.0_rk - poro(k), tiny)
                C_phase(k) = C_bulk(k) / one_m_phi
            end do
        end if
    end subroutine bulk_to_phase

    !------------------------------------------------------------
    ! Convert phase-specific -> bulk-sediment for ALL tracers
    ! Arrays are (nsed, ntr)
    !------------------------------------------------------------
    pure subroutine phase_to_bulk_all(C_phase, poro, tracer_info, C_bulk)
        real(rk), intent(in)  :: C_phase(:,:)          ! (nsed, ntr)
        real(rk), intent(in)  :: poro(:)               ! (nsed)
        type(TracerProperties), intent(in) :: tracer_info(:) ! (ntr)
        real(rk), intent(out) :: C_bulk(:,:)           ! (nsed, ntr)

        integer :: nsed, ntr, k, j
        real(rk) :: phi, one_m_phi
        real(rk), parameter :: tiny = 1.0e-12_rk

        nsed = size(C_phase, 1)
        ntr  = size(C_phase, 2)

        do j = 1, ntr
            if (tracer_info(j)%is_solute) then
                do k = 1, nsed
                    phi = max(poro(k), tiny)
                    C_bulk(k, j) = phi * C_phase(k, j)
                end do
            else
                do k = 1, nsed
                    one_m_phi = max(1.0_rk - poro(k), tiny)
                    C_bulk(k, j) = one_m_phi * C_phase(k, j)
                end do
            end if
        end do
    end subroutine phase_to_bulk_all

    !------------------------------------------------------------
    ! Convert bulk-sediment -> phase-specific for ALL tracers
    ! Arrays are (nsed, ntr)
    !------------------------------------------------------------
    pure subroutine bulk_to_phase_all(C_bulk, poro, tracer_info, C_phase)
        real(rk), intent(in)  :: C_bulk(:,:)           ! (nsed, ntr)
        real(rk), intent(in)  :: poro(:)               ! (nsed)
        type(TracerProperties), intent(in) :: tracer_info(:) ! (ntr)
        real(rk), intent(out) :: C_phase(:,:)          ! (nsed, ntr)

        integer :: nsed, ntr, k, j
        real(rk) :: phi, one_m_phi
        real(rk), parameter :: tiny = 1.0e-12_rk

        nsed = size(C_bulk, 1)
        ntr  = size(C_bulk, 2)

        do j = 1, ntr
            if (tracer_info(j)%is_solute) then
                do k = 1, nsed
                    phi = max(poro(k), tiny)
                    C_phase(k, j) = C_bulk(k, j) / phi
                end do
            else
                do k = 1, nsed
                    one_m_phi = max(1.0_rk - poro(k), tiny)
                    C_phase(k, j) = C_bulk(k, j) / one_m_phi
                end do
            end if
        end do
    end subroutine bulk_to_phase_all

    !! Clear sediment environment state and release memory.
    !! Deallocates all sediment arrays, clears the tridiagonal workspace, resets flags/counters,
    !! and nullifies the grid pointer.
    subroutine clear_sediment_env(SE)
        type(SedimentEnv), intent(inout) :: SE

        ! ---- Deallocate allocatable arrays ----
        if (allocated(SE%poro))       deallocate(SE%poro)
        if (allocated(SE%porewat_thickness)) deallocate(SE%porewat_thickness)
        if (allocated(SE%solid_thickness))  deallocate(SE%solid_thickness)
        if (allocated(SE%bioirr))     deallocate(SE%bioirr)

        if (allocated(SE%poro_w))   deallocate(SE%poro_w)
        if (allocated(SE%theta2))   deallocate(SE%theta2)
        if (allocated(SE%bioturb))  deallocate(SE%bioturb)
        if (allocated(SE%bioirr_w)) deallocate(SE%bioirr_w)

        if (allocated(SE%vel_solids))  deallocate(SE%vel_solids)
        if (allocated(SE%vel_solutes)) deallocate(SE%vel_solutes)

        if(allocated(SE%bulk_conc)) deallocate(SE%bulk_conc)
        if(allocated(SE%diff_sed)) deallocate(SE%diff_sed)
        if(allocated(SE%diff_sed0)) deallocate(SE%diff_sed0)
        if(allocated(SE%diff_sed_max)) deallocate(SE%diff_sed_max)
        if(allocated(SE%Db_eff_solids)) deallocate(SE%Db_eff_solids)
        if(allocated(SE%swi_flux)) deallocate(SE%swi_flux)

        ! ---- Reset workspace ----
        call clear_tridiag(SE%sed_trid)

        ! ---- Reset size/flags and pointer ----
        SE%nz      = 0
        SE%is_init = .false.
        nullify(SE%grid)
    end subroutine clear_sediment_env

    subroutine output_sed_profiles(full_grid, nsed, porosity_full, SE, static_prof)
        type(VerticalGrid),  intent(in)    :: full_grid
        integer,             intent(in)    :: nsed
        real(rk),            intent(in)    :: porosity_full(:)
        type(SedimentEnv),   intent(in)    :: SE
        type(StaticProfile), allocatable, intent(inout) :: static_prof(:)

        type(StaticProfile) :: prof
        integer :: nz_full
        real(rk), allocatable :: centre_full(:)
        real(rk), allocatable :: interface_full(:)

        nz_full = full_grid%nz

        if (nsed < 0 .or. nsed > nz_full) then
            error stop 'output_sed_profiles: invalid nsed.'
        end if

        if (size(porosity_full) /= nz_full) then
            error stop 'output_sed_profiles: porosity_full has wrong size.'
        end if

        if (nsed > 0) then
            if (.not. allocated(SE%theta2))      error stop 'output_sed_profiles: SE%theta2 not allocated.'
            if (.not. allocated(SE%vel_solids))  error stop 'output_sed_profiles: SE%vel_solids not allocated.'
            if (.not. allocated(SE%vel_solutes)) error stop 'output_sed_profiles: SE%vel_solutes not allocated.'
            if (SE%use_bioirrigation .and. .not. SE%output_bioirr_dynamic) then
                if (.not. allocated(SE%bioirr))  error stop 'output_sed_profiles: SE%bioirr not allocated.'
            end if
            if (SE%use_bioturbation .and. .not. SE%output_bioturb_dynamic) then
                if (.not. allocated(SE%bioturb)) error stop 'output_sed_profiles: SE%bioturb not allocated.'
            end if
        end if

        !============================================================
        ! Sediment mask on full grid (centres)
        !============================================================
        allocate(centre_full(nz_full))
        centre_full = 0.0_rk
        if (nsed > 0) centre_full(1:nsed) = 1.0_rk

        call init_static_profile(p            = prof, &
                                 name         = 'sediment_mask', &
                                 long_name    = 'Mask for sediment layers (1=sediment, 0=water)', &
                                 units        = '1', &
                                 profile_data = centre_full, &
                                 vert_coord   = 'centre' )
        prof%has_min   = .true.
        prof%valid_min = 0.0_rk
        prof%has_max   = .true.
        prof%valid_max = 1.0_rk
        call append_static_profile(static_prof, prof)
        deallocate(centre_full)

        !============================================================
        ! Porosity on full grid (centres)
        ! Water column naturally has porosity = 1
        !============================================================
        call init_static_profile(p            = prof, &
                                 name         = 'porosity', &
                                 long_name    = 'Porosity', &
                                 units        = '1', &
                                 profile_data = porosity_full, &
                                 vert_coord   = 'centre' )
        prof%has_min   = .true.
        prof%valid_min = 0.0_rk
        prof%has_max   = .true.
        prof%valid_max = 1.0_rk
        call append_static_profile(static_prof, prof)

        !============================================================
        ! Tortuosity squared on full grid (centres)
        ! Neutral extension into water column: 1
        !============================================================
        allocate(centre_full(nz_full))
        centre_full = 1.0_rk
        if (nsed > 0) centre_full(1:nsed) = SE%theta2(1:nsed)

        call init_static_profile(p            = prof, &
                                 name         = 'tortuosity2', &
                                 long_name    = 'Tortuosity squared', &
                                 units        = '1', &
                                 profile_data = centre_full, &
                                 vert_coord   = 'centre' )
        prof%has_min   = .true.
        prof%valid_min = 1.0_rk
        call append_static_profile(static_prof, prof)
        deallocate(centre_full)

        !============================================================
        ! Bioirrigation on full grid (centres)
        ! Zero in water column
        !============================================================
        if (SE%use_bioirrigation .and. .not. SE%output_bioirr_dynamic) then
            allocate(centre_full(nz_full))
            centre_full = 0.0_rk
            if (nsed > 0) centre_full(1:nsed) = SE%bioirr(1:nsed)

            call init_static_profile(p          = prof,            &
                                     name       = 'bioirrigation', &
                                     long_name  = 'Bioirrigation coefficient', &
                                     units      = 's-1', &
                                     profile_data = centre_full, &
                                     vert_coord = 'centre' )
            prof%has_min   = .true.
            prof%valid_min = 0.0_rk
            call append_static_profile(static_prof, prof)
            deallocate(centre_full)
        end if

        !============================================================
        ! Bioturbation on full grid (interfaces)
        ! Zero in water column
        !
        ! StaticProfile stores interface data as a packed 1..nz+1 array.
        ! Mapping from SE arrays with 0:nsed indexing:
        !   packed(i+1) <-> interface i
        !============================================================
        if (SE%use_bioturbation .and. .not. SE%output_bioturb_dynamic) then
            allocate(interface_full(nz_full + 1))
            interface_full = 0.0_rk
            if (nsed > 0) interface_full(1:nsed+1) = SE%bioturb(0:nsed)

            call init_static_profile(p          = prof, &
                                     name       = 'bioturbation', &
                                     long_name  = 'Bioturbation coefficient', &
                                     units      = 'm2 s-1', &
                                     profile_data = interface_full, &
                                     vert_coord = 'interface' )
            prof%has_min   = .true.
            prof%valid_min = 0.0_rk
            call append_static_profile(static_prof, prof)
            deallocate(interface_full)
        end if

        !============================================================
        ! Burial velocity of solids on full grid (interfaces)
        ! Zero in water column
        !============================================================
        allocate(interface_full(nz_full + 1))
        interface_full = 0.0_rk
        if (nsed > 0) interface_full(1:nsed+1) = SE%vel_solids(0:nsed)

        call init_static_profile(p          = prof, &
                                 name       = 'burial_velocity_solids', &
                                 long_name  = 'Burial velocity of solids', &
                                 units      = 'm s-1', &
                                 profile_data = interface_full, &
                                 vert_coord = 'interface' )
        call append_static_profile(static_prof, prof)
        deallocate(interface_full)

        !============================================================
        ! Burial velocity of solutes / porewater on full grid (interfaces)
        ! Zero in water column
        !============================================================
        allocate(interface_full(nz_full + 1))
        interface_full = 0.0_rk
        if (nsed > 0) interface_full(1:nsed+1) = SE%vel_solutes(0:nsed)

        call init_static_profile(p          = prof, &
                                 name       = 'burial_velocity_solutes', &
                                 long_name  = 'Burial velocity of porewater', &
                                 units      = 'm s-1', &
                                 profile_data = interface_full, &
                                 vert_coord = 'interface' )
        call append_static_profile(static_prof, prof)
        deallocate(interface_full)

    end subroutine output_sed_profiles



    !==================== SEDIMENT PROFILES=========================================
    !! Compute porosity and tortuosity profiles from an exponential compaction law.
    !! Returns porosity at layers centres (1:nz) and interfaces (0:nz), and squared
    !! tortuosity at interfaces using Boudreau (1997) and Soetaert (1997) formulation.
    subroutine compute_porosity_profile(grid, SedP, poro, poro_int, theta2, porewat_thickness, solid_thickness)
        type(VerticalGrid), intent(in) :: grid          ! this is the SEDIMENT grid
        type(SedParams),   intent(in)  :: SedP
        real(rk),          intent(out) :: poro(1:grid%nz)        ! Porosity at layers centres
        real(rk),          intent(out) :: poro_int(0:grid%nz)    ! Porosity at interfaces
        real(rk),          intent(out) :: theta2(0:grid%nz)      ! Diffusion tortuosity factor (Boudreau 1997), used as D_eff = D0 / theta2
        real(rk),          intent(out) :: porewat_thickness(1:grid%nz)  ! Storage capacity of porewater per unit horizontal-area
        real(rk),          intent(out) :: solid_thickness(1:grid%nz)   ! Storage capacity of solids per unit horizontal-area

        integer ::  k
        real(rk) :: z, por_decay 

        real(rk), parameter :: phi_min = 1.0e-6_rk
        real(rk), parameter :: phi_max = 1.0_rk - 1.0e-6_rk

        por_decay = max(SedP%poro_decay, 1.0e-12_rk)

        do k = 1, grid%nz   
            z = grid%z(k)   ! depth below SWI [m]
            poro(k) = SedP%poro_deep + (SedP%poro_sfc - SedP%poro_deep) * exp(-z / por_decay)
            poro(k) = min(max(poro(k), phi_min), phi_max)

            porewat_thickness(k) = poro(k) * grid%dz(k) 
            solid_thickness(k) =  (1.0_rk - poro(k)) * grid%dz(k) 
        end do

        do k = 0, grid%nz   
            z = grid%z_w(k) ! interfaces: z_w(nz)=0 at SWI
            ! Porosity at layer interfaces
            poro_int(k) = SedP%poro_deep + (SedP%poro_sfc - SedP%poro_deep) * exp(-z / por_decay)
            poro_int(k) = min(max(poro_int(k), phi_min), phi_max)  ! Ensuring it's within (0,1)

            !---Tortuosity (Boudreau 1997, Eq. 4.120) ------
            theta2(k) = 1.0_rk - 2.0_rk*log(poro_int(k))
        end do
    end subroutine compute_porosity_profile

    !! Compute bioturbation diffusivity profile at sediment interfaces.
    !! Uses a constant mixed layer (Db=db_surface) down to bioturbation_depth,
    !! followed by exponential decay with scale bioturbation_decay_depth.
    subroutine compute_bioturbation_profile(grid, db_surface, bioturbation_depth, bioturbation_decay_depth, Db_int)
        type(VerticalGrid), intent(in)  :: grid                     ! Sediment grid
        real(rk),           intent(in)  :: db_surface               ! [m2/s]
        real(rk),           intent(in)  :: bioturbation_depth       ! [m] thickness of well-mixed layer
        real(rk),           intent(in)  :: bioturbation_decay_depth ! [m] e-folding scale below that
        real(rk),           intent(out) :: Db_int(0:grid%nz)        ! [m2/s] interfaces (0:nz)

        integer :: nz, i
        real(rk) :: z, zdec

        nz = grid%nz        

        zdec = max(bioturbation_decay_depth, 1.0e-12_rk)

        ! Bioturbation profile according to Soetaert et al. (1996)
        do i = 0, nz
            z = grid%z_w(i)   ! [m] depth below SWI for sediment grid; SWI at z=0 is i=nz
            if (z <= bioturbation_depth) then
                Db_int(i) = db_surface
            else
                Db_int(i) = db_surface * exp(-(z - bioturbation_depth) / zdec)
            end if
        end do
    end subroutine compute_bioturbation_profile

    !! Compute bioirrigation exchange rate profile (alpha) at layers centres.
    !! Applies an exponential decay with depth below SWI, with surface rate alpha0
    !! and e-folding depth z_decay; outputs alpha(1:nz) and alpha_w(0:nz).
    !!
    !! Bioirrigation is represented as a non-local exchange term following Boudreau (1984),
    !! which arises from lateral averaging of burrow-scale irrigation processes. 
    !! The depth dependence of the exchange coefficient is parameterized as an exponential decay, 
    !! a commonly used empirical form motivated by observed decreases in burrow density
    !! and ventilation intensity with depth.
    subroutine compute_bioirrigation_alpha(grid, alpha0, z_decay, alpha, alpha_w)
        type(VerticalGrid), intent(in)  :: grid
        real(rk),           intent(in)  :: alpha0          ! [s-1]
        real(rk),           intent(in)  :: z_decay         ! [m]
        real(rk),           intent(out) :: alpha(1:grid%nz)  ! [s-1] (1:nz)
        real(rk),           intent(out) :: alpha_w(0:grid%nz)  ! [s-1] (0:nz)

        integer :: k, nz
        real(rk) :: z, zdec

        zdec = max(z_decay, 1.0e-12_rk)

        nz = grid%nz

        do k = 1, nz
            z = grid%z(k)    ! [m], positive downward. sediment-water interface at 0m
            alpha(k) = alpha0 * exp(-z/ zdec)
        end do

        ! Bioirrigation at layer the interfaces
        do k = 0, nz
            z = grid%z_w(k)    ! [m], positive downward. sediment-water interface at 0m
            alpha_w(k) = alpha0 * exp(-z/ zdec)
        end do
    end subroutine compute_bioirrigation_alpha

    !===================================================================
    !         Burial velocities for solids and pore water 
    !     under the assumption of steady-state compaction
    ! Assumes steady-state compaction with constant solid volume flux
    ! (1-phi)w = (1-phi_inf)w_inf and constant porewater volume flux
    ! phi*u = phi_inf*w_inf, where phi_inf is the asymptotic porosity
    ! and w_inf is the deep solid burial velocity.
    !==================================================================
    subroutine compute_burial_velocities(grid, SedP, poro_int, w_inf, w_solid, u_pore)
        type(VerticalGrid), intent(in)  :: grid
        type(SedParams),    intent(in)  :: SedP
        real(rk), intent(in)  :: poro_int(0:grid%nz)   ! porosity at interfaces
        real(rk), intent(in)  :: w_inf                 ! Solid burial velocity at infinite depth. [m/s]
        real(rk), intent(out) :: w_solid(0:grid%nz)    ! Solid-phase burial velocity [m/s]
        real(rk), intent(out) :: u_pore(0:grid%nz)     ! Porewater burial velocity [m/s]

        integer  :: k
        real(rk) :: phi_inf, phik
        real(rk) :: winf_neg
        real(rk), parameter :: phi_min = 1.0e-12_rk

        winf_neg = - abs(w_inf)      ! Enforcing negative velocity for burial        

        ! Asymptotic porosity (CANDI-style reference)
        phi_inf = min(max(SedP%poro_deep, phi_min), 1.0_rk - phi_min) 

        do k = 0, grid%nz
            phik = min(max(poro_int(k), phi_min), 1._rk - phi_min)

            ! Solid burial velocity: (1-phi) w  = (1-phi_inf) winf_neg
            ! (1-phi) w = (1-phi_inf) w_inf   (Boudreau, 1997, Eq. 3.67)
            w_solid(k) = winf_neg * (1._rk - phi_inf) / (1._rk - phik)

            ! Porewater burial velocity: phi u = phi_inf winf_neg
            ! phi u = phi_inf w_inf           (Boudreau, 1997, Eq. 3.68)
            u_pore(k)  = (winf_neg * phi_inf) / phik
        end do
        ! Enforce no advective burial flux through the SWI (top interface k=nz).
        ! SWI exchange is handled separately via deposition/solute exchange.
        !w_solid(nz) = 0.0_rk
        !u_pore(nz) = 0.0_rk
    end subroutine compute_burial_velocities

    !================= Helpers================== 

    !! Convert sediment parameters from units commonly reported in the literature to SI units.
    !! Depths: cm->m; rates: cm/yr->m/s, cm2/yr->m2/s, 1/yr->1/s. Porosities are unchanged.
    pure subroutine convert_units_to_SI(user, si)
        type(SedParams), intent(in)  :: user   ! user units (cm, yr)
        type(SedParams), intent(out) :: si     ! SI units (m, s)

        real(rk), parameter :: cm_to_m = 1.0e-2_rk
        real(rk), parameter :: yr_to_s = 365.25_rk * 86400.0_rk

        ! Copy porosity params (dimensionless)
        si%poro_sfc  = user%poro_sfc
        si%poro_deep = user%poro_deep
        si%biot_mode = to_lower(user%biot_mode)
        si%irr_mode  = to_lower(user%irr_mode)

        ! Depth scales: cm -> m
        si%poro_decay    = user%poro_decay    * cm_to_m
        si%biot_mld      = user%biot_mld      * cm_to_m
        si%biot_ez       = user%biot_ez       * cm_to_m
        si%irr_ez        = user%irr_ez        * cm_to_m

        ! Rates:
        ! sed_rate: cm/yr -> m/s
        si%sed_rate = user%sed_rate * cm_to_m / yr_to_s

        ! Bioturbation coefficient: cm^2/yr -> m^2/s
        si%biot_db_sfc = user%biot_db_sfc * (cm_to_m*cm_to_m) / yr_to_s

        ! irr_sfc: 1/yr -> 1/s
        si%irr_sfc = user%irr_sfc / yr_to_s
        ! Copy cnpar value
        si%cnpar_sed = user%cnpar_sed
    end subroutine convert_units_to_SI

    ! Allocate arrays to store sediment properties
    subroutine allocate_sed_arrays(SE)
        type(SedimentEnv), intent(inout) :: SE
        integer :: nz

        if (.not. associated(SE%grid)) error stop "allocate_sed_arrays: SE%grid not associated"
        nz = SE%grid%nz

        ! Centres: 1:nz        
        if (allocated(SE%poro))   deallocate(SE%poro)
        allocate(SE%poro(nz))

        if (allocated(SE%porewat_thickness))   deallocate(SE%porewat_thickness)
        allocate(SE%porewat_thickness(nz))

        if (allocated(SE%solid_thickness))   deallocate(SE%solid_thickness)
        allocate(SE%solid_thickness(nz))
               
        if (allocated(SE%bioirr)) deallocate(SE%bioirr)
        allocate(SE%bioirr(nz))        

        ! Interfaces: 0:nz  -> size = nz+1        
        if (allocated(SE%poro_w)) deallocate(SE%poro_w)
        allocate(SE%poro_w(0:nz))
        
        if (allocated(SE%theta2)) deallocate(SE%theta2)
        allocate(SE%theta2(0:nz))

        if (allocated(SE%bioturb)) deallocate(SE%bioturb)
        allocate(SE%bioturb(0:nz))

        if (allocated(SE%bioirr_w)) deallocate(SE%bioirr_w)
        allocate(SE%bioirr_w(0:nz)) 

        ! Array to store molecular diffusivities in the sediments (unique per tracer)
        if (allocated(SE%diff_sed)) deallocate(SE%diff_sed)
        allocate(SE%diff_sed(0:nz)) 
        SE%diff_sed = 0.0_rk

        ! Array to store maximum molecular diffusivities in the sediments 
        if (allocated(SE%diff_sed_max)) deallocate(SE%diff_sed_max)
        allocate(SE%diff_sed_max(0:nz)) 
        SE%diff_sed_max = 0.0_rk


        if (allocated(SE%Db_eff_solids)) deallocate(SE%Db_eff_solids)
        allocate(SE%Db_eff_solids(0:nz)) 
        SE%Db_eff_solids = 0.0_rk

        ! --- Burial velocities (interfaces 0:nz)
        if (allocated(SE%vel_solids)) deallocate(SE%vel_solids)
        allocate(SE%vel_solids(0:nz))
        SE%vel_solids = 0._rk

        if (allocated(SE%vel_solutes)) deallocate(SE%vel_solutes)
        allocate(SE%vel_solutes(0:nz))
        SE%vel_solutes = 0._rk
    end subroutine allocate_sed_arrays

    !! Write sediment interface profiles to an ASCII diagnostics file.
    !! Outputs a single table on the interface grid (0:nz): depth, porosity, tortuosity,
    !! bioturbation diffusivity, bioirrigation exchange rate, and burial velocities.
    subroutine write_sediment_profiles(filename, grid, SE, ierr)
        character(*),       intent(in)  :: filename
        type(VerticalGrid), intent(in)  :: grid
        type(SedimentEnv),  intent(in)  :: SE
        integer, optional,  intent(out) :: ierr

        integer :: u, k, nz

        character(len=*), parameter :: fmt_row = &
            '(I6,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10)'

        if (present(ierr)) ierr = 0
        nz = grid%nz

        if (.not. allocated(SE%poro_w) .or. .not. allocated(SE%theta2) .or. &
            .not. allocated(SE%bioturb) .or. .not. allocated(SE%bioirr_w) .or. &
            .not. allocated(SE%vel_solids) .or. .not. allocated(SE%vel_solutes)) then
            if (present(ierr)) then
                ierr = 2
                return
            else
                error stop "write_sediment_profiles: Missing sediment interface arrays"
            end if
        end if

        open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

        write(u,'(A)') '# Sediment interface profiles'
        write(u,'(A)') '# z_w is depth below SWI [m], positive downward; z_w(nz)=0 at SWI; z_w(0)=sediment bottom'
        write(u,'(A)') '# Velocities are burial velocities at layer interfaces; negative values indicate downward motion.'
        write(u,'(A)') 'Layer   depth[m]  porosity[-]  tortuosity[squared]  Db[m2s-1]  bioirr[s-1]  w_solid[ms-1]  u_pore[ms-1]'

        ! Print from SWI downward for readability (k = nz .. 0)
        do k = nz, 0, -1
            write(u,fmt_row) k, grid%z_w(k), SE%poro_w(k), SE%theta2(k), SE%bioturb(k), SE%bioirr_w(k), &
                   SE%vel_solids(k), SE%vel_solutes(k)
        end do

        close(u)
    end subroutine write_sediment_profiles


    ! Writes a file with the tracer properties
    subroutine write_tracer_properties(BE, filename)
        type(BioEnv),     intent(in) :: BE
        character(len=*), intent(in) :: filename

        integer :: iu, i, ntr
        character(len=64) :: name
        character(len=11) :: phase
        character(len=24) :: dm

        ! Strings for NA / numeric columns
        character(len=16) :: s_ads
        character(len=16) :: sA, sB, sm0, sm1, sA0, sEa, sVb, sDref, sTref

        if (.not. allocated(BE%tracer_info)) then
            write(*,*) 'write_tracer_properties: tracer_info not allocated.'
            return
        end if

        ntr = size(BE%tracer_info)
        if (ntr <= 0) return

        open(newunit=iu, file=trim(filename), status='replace', action='write', form='formatted')

        write(iu,'(a)') '# idx name phase adsorp diff_method A B m0 m1 A0 Ea Vb Dref Tref'
        write(iu,'(a)') '# NA = not applicable'
        write(iu,'(a)') 'idx  name           phase        adsorp           diff_method                 ' // &
                        'A                B               m0               m1               ' // &
                        'A0               Ea               Vb              Dref             Tref'
        write(iu,'(a)') '---  -------------  -----------  --------------  -------------------------  ' // &
                        '---------------  ---------------  ---------------  ---------------  ' // &
                        '---------------  ---------------  ---------------  ---------------  ---------------'


        do i = 1, ntr
            ! Name
            name = 'unknown'
            if (associated(BE%model)) then
                if (size(BE%model%interior_state_variables) >= i) then
                    name = trim(BE%model%interior_state_variables(i)%name)
                end if
            end if

            ! Phase
            if (BE%tracer_info(i)%is_solute) then
                phase = 'solute'
            else
                phase = 'particulate'
            end if

            ! Adsorption (NA for particulates)
            if (BE%tracer_info(i)%is_solute) then
                s_ads = fmt_real(BE%tracer_info(i)%adsorption)
            else
                s_ads = 'NA'
            end if

            ! Diff method label
            dm = diff_label(BE%tracer_info(i)%diff_method)

            ! Default all params to NA
            sA    = 'NA'; sB    = 'NA'
            sm0   = 'NA'; sm1   = 'NA'
            sA0   = 'NA'; sEa   = 'NA'
            sVb   = 'NA'
            sDref = 'NA'; sTref = 'NA'

            ! Populate only the relevant parameter columns
            select case (BE%tracer_info(i)%diff_method)
            case (DIFF_O2CO2_AB)
                sA = fmt_real(BE%tracer_info(i)%A)
                sB = fmt_real(BE%tracer_info(i)%B)

            case (DIFF_ION_LINEAR)
                sm0 = fmt_real(BE%tracer_info(i)%m0)
                sm1 = fmt_real(BE%tracer_info(i)%m1)

            case (DIFF_ARRHENIUS)
                sA0 = fmt_real(BE%tracer_info(i)%A0)
                sEa = fmt_real(BE%tracer_info(i)%Ea)

            case (DIFF_WILKE_CHANG)
                sVb = fmt_real(BE%tracer_info(i)%Vb)

            case (DIFF_STOKES_EINSTEIN)
                sDref = fmt_real(BE%tracer_info(i)%Dref)
                sTref = fmt_real(BE%tracer_info(i)%Tref)

            case default
                ! DIFF_NONE or unknown -> leave as NA
            end select

            ! Write row (IMPORTANT: widen a-fields so values never truncate)
            write(iu,'(i3,2x,a13,2x,a11,2x,a14,2x,a25, 2x,a15,2x,a15,2x,a15,2x,a15,2x,a15,2x,a15,2x,a15,2x,a15,2x,a15)') &
                BE%tracer_info(i)%fabm_index, trim(name), trim(phase), trim(s_ads), trim(dm), &
                trim(sA), trim(sB), trim(sm0), trim(sm1), trim(sA0), trim(sEa), trim(sVb), trim(sDref), trim(sTref)
        end do

        close(iu)

    contains

        pure function diff_label(method) result(s)
            integer, intent(in) :: method
            character(len=24) :: s
            select case (method)
            case (DIFF_NONE);            s = 'DIFF_NONE(0)'
            case (DIFF_O2CO2_AB);        write(s,'(a,i0,a)') 'DIFF_O2CO2_AB(', method, ')'
            case (DIFF_ION_LINEAR);      write(s,'(a,i0,a)') 'DIFF_ION_LINEAR(', method, ')'
            case (DIFF_ARRHENIUS);       write(s,'(a,i0,a)') 'DIFF_ARRHENIUS(', method, ')'
            case (DIFF_WILKE_CHANG);     write(s,'(a,i0,a)') 'DIFF_WILKE_CHANG(', method, ')'
            case (DIFF_STOKES_EINSTEIN); write(s,'(a,i0,a)') 'DIFF_STOKES_EINSTEIN(', method, ')'
            case default;                write(s,'(a,i0,a)') 'DIFF_UNKNOWN(', method, ')'
            end select
        end function diff_label

        pure function fmt_real(x) result(s)
            real(rk), intent(in) :: x
            character(len=16) :: s
            real(rk) :: ax
            ax = abs(x)

            if (ax == 0._rk) then
                s = '0'
            else if (ax >= 1.e-3_rk .and. ax < 1.e4_rk) then
                ! "Normal-looking" decimals; 6 dp is a decent default for these tables
                write(s,'(G0.6)') x
            else
                ! Tiny/huge -> scientific
                write(s,'(ES12.4)') x
            end if

            s = adjustl(s)
        end function fmt_real

    end subroutine write_tracer_properties    

end module sediment