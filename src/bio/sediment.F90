module sediment    
    use bio_types,         only: SedimentEnv
    use bio_params,        only: SedParams, read_sed_parameters
    use grids,             only: VerticalGrid
    use precision_types,   only: rk
    use read_config_yaml,  only: ConfigParams
    use tridiagonal,       only: init_tridiag, clear_tridiag

    implicit none
    private

    public :: init_sediment, clear_sediment_env

contains

    ! Initialise sediment environment for one column.
    !! Reads user parameters from YAML, converts to SI units, allocates arrays,
    !! and precomputes sediment property profiles (porosity, tortuosity, Db, irrigation).
    subroutine init_sediment(cfg_params, grid, SE)
        type(ConfigParams),  intent(in) :: cfg_params
        type(VerticalGrid), intent(in), target  :: grid        ! Sediment grid
        type(SedimentEnv),intent(inout) :: SE

        ! Associate grid pointer
        SE%grid => grid
        SE%nz = grid%nz

        ! Read sediment parameters
        call read_sed_parameters(cfg_params, grid, SE%params_user)
        ! Convert parameter units to SI units cm->m and yr->s
        call convert_units_to_SI(SE%params_user, SE%params_SI)
        ! Allocate working arrays 
        call allocate_sed_arrays(SE)

        ! Compute profiles for sediment properties
        call compute_porosity_profile(grid, SE%params_SI, SE%poro, SE%poro_w, SE%theta)
        call compute_bioturbation_profile(grid, SE%params_SI%biot_db_sfc, SE%params_SI%biot_mld, SE%params_SI%biot_ez, SE%bioturb)
        call compute_bioirrigation_alpha(grid, SE%params_SI%irr_sfc, SE%params_SI%irr_ez, SE%bioirr, SE%bioirr_w)

!!! Initialise tridiagonal matrix (Later)

        SE%is_init = .true.
        call write_sediment_profiles('Sediment_profiles.dat',SE%grid, SE)
    end subroutine init_sediment

    !! Compute porosity and tortuosity profiles from an exponential compaction law.
    !! Returns porosity at layers centres (1:nz) and interfaces (0:nz), and squared
    !! tortuosity at interfaces using Boudreau (1997) and Soetaert (1997) formulation.
    subroutine compute_porosity_profile(grid, SedP, poro, poro_int, theta2)
        type(VerticalGrid), intent(in) :: grid          ! this is the SEDIMENT grid
        type(SedParams),   intent(in)  :: SedP
        real(rk),          intent(out) :: poro(1:grid%nz)      ! Porosity at layers centres
        real(rk),          intent(out) :: poro_int(0:grid%nz)  ! Porosity at interfaces
        real(rk),          intent(out) :: theta2(0:grid%nz)    ! Squared tortuosity at interfaces

        integer ::  k
        real(rk) :: z, por_decay, phi   

        por_decay = max(SedP%poro_decay, 1.0e-12_rk)

        do k = 1, grid%nz   
            z = grid%z(k)   ! depth below SWI [m]
            poro(k) = SedP%poro_deep + (SedP%poro_sfc - SedP%poro_deep) * exp(-z / por_decay)
        end do

        do k = 0, grid%nz   
            z = grid%z_w(k) ! interfaces: z_w(nz)=0 at SWI
            ! Porosity at layer interfaces
            poro_int(k) = SedP%poro_deep + (SedP%poro_sfc - SedP%poro_deep) * exp(-z / por_decay)

            phi = min(max(poro_int(k), 1.0e-6_rk), 1.0_rk - 1.0e-6_rk)  ! Ensuring it's within (0,1)
            !---Tortuosity (Boudreau 1997, Eq. 4.120) ------
            theta2(k) = 1.0_rk - 2.0_rk*log(phi)
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
        real(rk),           intent(out) :: Db_int(0:grid%nz)             ! [m2/s] interfaces (0:nz)

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
    subroutine compute_bioirrigation_alpha(grid, alpha0, z_decay, alpha, alpha_w)
        type(VerticalGrid), intent(in)  :: grid
        real(rk),           intent(in)  :: alpha0          ! [s-1]
        real(rk),           intent(in)  :: z_decay         ! [m]
        real(rk),           intent(out) :: alpha(grid%nz)  ! [s-1] (1:nz)
        real(rk),           intent(out) :: alpha_w(0:grid%nz)  ! [s-1] (1:nz)

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
    !! Applies constant solid flux (1-phi)w and porewater flux phi*u relative to a
    !! reference depth z_ref and burial velocity w_ref (Boudreau, 1997).
    !==================================================================
    subroutine compute_burial_velocities(grid, poro_int, z_ref, w_ref, w_solid, u_pore)
        type(VerticalGrid), intent(in)  :: grid
        real(rk), intent(in)  :: poro_int(0:grid%nz)   ! porosity at interfaces
        real(rk), intent(in)  :: z_ref                 ! [m] reference depth below SWI
        real(rk), intent(in)  :: w_ref                 ! [m/s] burial velocity at z_ref
        real(rk), intent(out) :: w_solid(0:grid%nz)    ! Burial velocity of particulate matter [m/s]
        real(rk), intent(out) :: u_pore(0:grid%nz)     ! Burial velocity of pore water [m/s]

        integer  :: k, nz, k_ref
        real(rk) :: phi_ref, phik, phi_min

        nz = grid%nz
        phi_min = 1.0e-12_rk

        ! Find index for closest interface to reference depth
        k_ref = minloc(abs(grid%z_w(0:nz) - z_ref), dim=1) - 1
        phi_ref = min(max(poro_int(k_ref), phi_min), 1._rk - phi_min)

        do k = 0, nz
            phik = min(max(poro_int(k), phi_min), 1._rk - phi_min)

            ! Solid burial velocity: (1-phi) w  = (1-phi_ref) w_ref
            ! (1-ϕ)w = (1-ϕx)wx  Eq. 3.67 in Boudreau (1997)
            w_solid(k) = w_ref * (1._rk - phi_ref) / (1._rk - phik)

            ! Porewater burial velocity: phi u = phi_ref w_ref
            ! ϕu = ϕx wx     Eq. 3.68 in Boudreau (1997)
            u_pore(k)  = (w_ref * phi_ref) / phik
        end do
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

        ! Depth scales: cm -> m
        si%sed_ref_depth = user%sed_ref_depth * cm_to_m
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
               
        if (allocated(SE%bioirr)) deallocate(SE%bioirr)
        allocate(SE%bioirr(nz))        

        ! Interfaces: 0:nz  -> size = nz+1        
        if (allocated(SE%poro_w)) deallocate(SE%poro_w)
        allocate(SE%poro_w(0:nz))
        
        if (allocated(SE%theta)) deallocate(SE%theta)
        allocate(SE%theta(0:nz))

        if (allocated(SE%bioturb)) deallocate(SE%bioturb)
        allocate(SE%bioturb(0:nz))

        if (allocated(SE%bioirr_w)) deallocate(SE%bioirr_w)
        allocate(SE%bioirr_w(0:nz)) 
    end subroutine allocate_sed_arrays

    !! Write sediment interface profiles to an ASCII diagnostics file.
    !! Outputs a single table on the interface grid (0:nz): depth, porosity, tortuosity^2,
    !! bioturbation diffusivity, and bioirrigation exchange rate.
    subroutine write_sediment_profiles(filename, grid, SE, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(*),       intent(in)  :: filename
        type(VerticalGrid), intent(in)  :: grid
        type(SedimentEnv),  intent(in)  :: SE
        integer, optional,  intent(out) :: ierr

        integer :: u, k, nz
        if (present(ierr)) ierr = 0

        nz = grid%nz

        if (.not. allocated(SE%poro_w) .or. .not. allocated(SE%theta) .or. &
            .not. allocated(SE%bioturb) .or. .not. allocated(SE%bioirr_w)) then
            if (present(ierr)) then
                ierr = 2; return
            else
                error stop "write_sediment_profiles: Missing sediment interface arrays"
            end if
        end if

        ! Size checks (avoid silent mismatch)
        if (size(SE%poro_w) /= nz+1 .or. size(SE%theta) /= nz+1 .or. &
            size(SE%bioturb) /= nz+1 .or. size(SE%bioirr_w) /= nz+1) then
            if (present(ierr)) then
                ierr = 3; return
            else
                error stop "write_sediment_profiles: Array size mismatch (expected 0:nz)"
            end if
        end if

        ! ---- Write file ----
        open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

        write(u,'(A)') '# Sediment interface profiles'
        write(u,'(A)') '# Depth coordinate: z_w is depth below SWI [m], positive downward; z_w(nz)=0 at SWI'
        write(u,'(A)') 'Layer   depth[m]  porosity[-]  tortuosity[squared]  Db[m2s-1]  bioirr[s-1]'

        ! Print from SWI downward for readability (k = nz .. 0)
        do k = nz, 0, -1
            write(u,'(I6,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10,1X,ES18.10)') &
                k, grid%z_w(k), SE%poro_w(k), SE%theta(k), SE%bioturb(k), SE%bioirr_w(k)
        end do

        close(u)
    end subroutine write_sediment_profiles


    !! Clear sediment environment state and release memory.
    !! Deallocates all sediment arrays, clears the tridiagonal workspace, resets flags/counters,
    !! and nullifies the grid pointer.
    subroutine clear_sediment_env(SE)
        type(SedimentEnv), intent(inout) :: SE

        ! ---- Deallocate allocatable arrays ----
        if (allocated(SE%poro))     deallocate(SE%poro)
        if (allocated(SE%bioirr))   deallocate(SE%bioirr)

        if (allocated(SE%poro_w))   deallocate(SE%poro_w)
        if (allocated(SE%theta))    deallocate(SE%theta)
        if (allocated(SE%bioturb))  deallocate(SE%bioturb)
        if (allocated(SE%bioirr_w)) deallocate(SE%bioirr_w)

        ! ---- Reset workspace ----
        call clear_tridiag(SE%trid)

        ! ---- Reset size/flags and pointer ----
        SE%nz      = 0
        SE%is_init = .false.
        nullify(SE%grid)
    end subroutine clear_sediment_env




end module sediment