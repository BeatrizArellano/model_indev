module bio_main
    use fabm
    use bio_params,          only: read_bio_parameters, mol_diff
    use bio_types,           only: BioEnv
    use geo_utils,           only: LocationInfo
    use physics_types,       only: PhysicsState
    use forcing_manager,     only: ForcingSnapshot
    use grids,               only: VerticalGrid   
    use light,               only: compute_SW_profile, compute_PAR_profile
    use numerical_stability, only: compute_bio_substeps
    use pressure,            only: compute_pressure 
    use precision_types,     only: rk, lk
    use read_config_yaml,    only: ConfigParams
    use tridiagonal,         only: init_tridiag, clear_tridiag
    use vertical_transport,  only: apply_vertical_transport
    use vertical_mixing,     only: scalar_diffusion


  implicit none
  private

  public :: init_bio_fabm, integrate_bio_fabm, end_bio_fabm


contains

    subroutine init_bio_fabm(cfg, location, grid, timestep, PS, FS, BE)
        type(ConfigParams),       intent(in) :: cfg
        type(LocationInfo),       intent(in) :: location
        type(VerticalGrid),       intent(in) :: grid
        integer(lk), intent(in)              :: timestep
        type(PhysicsState),       intent(in) :: PS
        type(ForcingSnapshot),    intent(in) :: FS
        type(BioEnv),          intent(inout) :: BE

        integer  :: maxlen, ivar, nz, nint, nsfc, nbtm
        real(rk) :: dt_main

        write(*,'(A)') 'Initialising biogeochemistry via FABM...'

        BE%grid     = grid   ! Full grid
        BE%wat_grid = grid   ! For now they're the same grid
        dt_main     = real(timestep, rk)
        nz          = BE%grid%nz              ! Number of vertical layers 

        call read_bio_parameters(cfg,BE%params)

        ! Initialize (reads FABM configuration from fabm.yaml)
        ! After this the number of biogeochemical variables is fixed.
        BE%model => fabm_create_model(BE%params%config_file)
    

        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        call BE%model%set_domain(nz, seconds_per_time_unit=dt_main)

        ! ---- Interior state variables ---------------
        nint = size(BE%model%interior_state_variables)
        BE%BS%n_interior = nint
        if (nint > 0) then
            ! Allocate the array for interior tracers concentrations
            if(allocated(BE%BS%interior_state)) deallocate(BE%BS%interior_state)
            allocate(BE%BS%interior_state(nz, nint))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%interior_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%intvar_names(nint))
            ! Point FABM to the current concentration state for biogeochemical tracers
            do ivar = 1, nint
                BE%BS%intvar_names(ivar) = trim(BE%model%interior_state_variables(ivar)%name) ! Retrieve names of tracers
                call BE%model%link_interior_state_data(ivar, BE%BS%interior_state(:,ivar))
            end do  
        end if
        ! Allocating working arrays        
        allocate(BE%tendency_int(nz,nint))                        ! Sources (dC/dt) 
        allocate(BE%velocity(nz,nint))                            ! Velocities
        allocate(BE%flux_sf(nint));    allocate(BE%flux_bt(nint)) ! Fluxes at the surface and bottom of interior variables 


        ! --- Bottom-only state variables ---
        nbtm = size(BE%model%bottom_state_variables)
        BE%BS%n_bottom = nbtm
        if (nbtm > 0) then
            if (allocated(BE%BS%bottom_state)) deallocate(BE%BS%bottom_state)
            allocate(BE%BS%bottom_state(nbtm))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%bottom_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%btmvar_names(nbtm))
            ! Link FABM to the bottom state vector
            call BE%model%link_all_bottom_state_data(BE%BS%bottom_state)
            do ivar=1, nbtm
                BE%BS%btmvar_names(ivar) = trim(BE%model%bottom_state_variables(ivar)%name) ! Retrieve names of tracers
            end do
        end if
        ! Allocate working array for sources (dC/dt) at the bottom
        allocate(BE%tendency_bt(nbtm))

        ! --- Surface-only state variables ---
        nsfc = size(BE%model%surface_state_variables)
        BE%BS%n_surface = nsfc
        if (nsfc > 0) then
            if (allocated(BE%BS%surface_state)) deallocate(BE%BS%surface_state)
            allocate(BE%BS%surface_state(nsfc))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%surface_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%sfcvar_names(nsfc))
            ! Link FABM to the surface state vector
            call BE%model%link_all_surface_state_data(BE%BS%surface_state)
            do ivar=1, nsfc
                BE%BS%sfcvar_names(ivar) = trim(BE%model%surface_state_variables(ivar)%name) ! Retrieve names of tracers
            end do
        end if
        ! Allocate working array for sources (dC/dt) at the surface
        allocate(BE%tendency_sf(nsfc));

                   
             
        
        BE%BS%n_total = nint + nsfc + nbtm  ! Total number of variables

        ! Provide location info to FABM
        call BE%model%link_interior_data(fabm_standard_variables%depth,BE%grid%z)                ! layer depths at centres
        call BE%model%link_interior_data(fabm_standard_variables%cell_thickness,BE%grid%dz)      ! layer thicknesses
        call BE%model%link_horizontal_data(fabm_standard_variables%latitude, location%lat)       ! Latitude
        call BE%model%link_horizontal_data(fabm_standard_variables%longitude, location%lon)      ! Longitude
        call BE%model%link_horizontal_data(fabm_standard_variables%bottom_depth, location%depth) ! Bottom depth
        ! call model%set_bottom_index(bottom_indices)
        ! call model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, yearday)

        ! Getting runtime variable ids to link environmental data to FABM.
        BE%id_temp    = BE%model%get_interior_variable_id(fabm_standard_variables%temperature)
        BE%id_salt    = BE%model%get_interior_variable_id(fabm_standard_variables%practical_salinity)
        BE%id_rho     = BE%model%get_interior_variable_id(fabm_standard_variables%density)
        BE%id_par     = BE%model%get_interior_variable_id(fabm_standard_variables%downwelling_photosynthetic_radiative_flux)
        BE%id_swr     = BE%model%get_interior_variable_id(fabm_standard_variables%downwelling_shortwave_flux)
        BE%id_pres    = BE%model%get_interior_variable_id(fabm_standard_variables%pressure)
        ! Horizontal variables
        BE%id_windspd = BE%model%get_horizontal_variable_id(fabm_standard_variables%wind_speed)
        BE%id_slp     = BE%model%get_horizontal_variable_id(fabm_standard_variables%surface_air_pressure)
        BE%id_par_sfc = BE%model%get_horizontal_variable_id(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        BE%id_swr_sfc = BE%model%get_horizontal_variable_id(fabm_standard_variables%surface_downwelling_shortwave_flux)
        BE%id_cloud   = BE%model%get_horizontal_variable_id(fabm_standard_variables%cloud_area_fraction)
        BE%id_stressb = BE%model%get_horizontal_variable_id(fabm_standard_variables%bottom_stress)
        BE%id_co2     = BE%model%get_horizontal_variable_id(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air)
        
        ! Finding out which interior variables are actually used in FABM
        BE%need_temp    = BE%model%is_variable_used(BE%id_temp) 
        BE%need_salt    = BE%model%is_variable_used(BE%id_salt)
        BE%need_rho     = BE%model%is_variable_used(BE%id_rho)
        BE%need_par     = BE%model%is_variable_used(BE%id_par)
        BE%need_swr     = BE%model%is_variable_used(BE%id_swr)
        BE%need_pres    = BE%model%is_variable_used(BE%id_pres) .or. BE%model%variable_needs_values(BE%id_pres)
        BE%need_windspd = BE%model%is_variable_used(BE%id_windspd)
        BE%need_slp     = BE%model%is_variable_used(BE%id_slp) 
        BE%need_par_sfc = BE%model%is_variable_used(BE%id_par_sfc)
        BE%need_swr_sfc = BE%model%is_variable_used(BE%id_swr_sfc)
        BE%need_cloud   = BE%model%is_variable_used(BE%id_cloud)
        BE%need_stressb = BE%model%is_variable_used(BE%id_stressb)
        BE%need_co2     = BE%model%is_variable_used(BE%id_co2) .or. BE%model%variable_needs_values(BE%id_co2)


        if (BE%need_pres) then
            if (.not. allocated(BE%BS%pres)) allocate(BE%BS%pres(nz))
        end if
        if (BE%need_par .and. .not. allocated(BE%BS%par)) then
            allocate(BE%BS%par(nz))
        end if
        if (BE%need_swr .and. .not. allocated(BE%BS%swr)) then
            allocate(BE%BS%swr(nz))
        end if
        if (BE%need_cloud) then
            write(*,*) 'Access to cloud cover data is not already implemented'
            stop 1
        end if

        ! Check other potential variables needed -> Implement a way to provide those variables. 
        ! Loop over a section within biogeochemistry in the yaml file to load those variables. 

    !!! Change when sediments are implemented   
        ! Point FABM to environmental data    
        call link_environment_data(PS, FS, BE)      

        ! Complete initialization and check whether FABM has all dependencies fulfilled
        ! (i.e., whether all required calls to model%link_*_data have been made)
        call BE%model%start()

    !!! Later implement initialising from a previous state loaded via a file. 
        ! Initialize the tracers
        ! This sets the values of the arrays: interior_state, bottom_state and surface_state
        if (nint > 0) call BE%model%initialize_interior_state(1,nz)
        if (nsfc > 0) call BE%model%initialize_surface_state()
        if (nbtm > 0) call BE%model%initialize_bottom_state()

        if (BE%BS%n_total <= 0) then
            write(*,*) "No variables found in biogeochemistry, disable biogeochemistry or provide a valid configuration of modules."
            stop 1
        end if

        ! Tridiagonal workspace
        call init_tridiag(BE%trid, nz)
        
        BE%is_init = .true.
        write(*,'("✓ Biogeochemistry initialised successfully with ",I0," variables:")') BE%BS%n_total
        write(*,'(" ",I0," interior, ",I0," surface and ",I0," bottom variables")') BE%BS%n_interior, BE%BS%n_surface, BE%BS%n_bottom
    end subroutine init_bio_fabm

    !=====================================================================================
    ! Update biogeochemistry for one timestep.
    !   - Communicates to FABM current environmental state.
    !   - Gets from FABM source terms and boundary fluxes.
    !   - Applies vertical mixing and moves tracers according to reported velocities.
    !   - Integrates tracer tendencies (sources: dC/dt) subcycling if needed for stability.
    !=====================================================================================
    subroutine integrate_bio_fabm(BE, PS, FS, timestep, istep_main)
        type(BioEnv),          intent(inout) :: BE
        type(PhysicsState),    intent(in)    :: PS
        type(ForcingSnapshot), intent(in)    :: FS
        integer(lk),           intent(in)    :: timestep
        integer(lk),           intent(in)    :: istep_main

        integer  :: nz, ivar, k, nint, nsfc, nbtm
        real(rk) :: istep_rk, dt_main, dt_sub
        integer  :: nsub, isub
        integer :: ierr          

        if (.not. BE%is_init) then
            write(*,*) 'integrate_bio_fabm: Biogeochemistry is not initialised.'
            return
        end if

        nz       = BE%grid%nz              ! Number of vertical layers       
        istep_rk = real(istep_main, rk)    ! Current time-step number in the main loop
        dt_main  = real(timestep, rk)
        nint     = BE%BS%n_interior        ! Number of interior tracers
        nsfc     = BE%BS%n_surface
        nbtm     = BE%BS%n_bottom;

        !---------------------------------------------
        ! Update FABM with environment data
        !--------------------------------------------
        call link_environment_data(PS, FS, BE)

        ! Repair state, if needed, to let FABM restore all state variables within their registered bounds.
        call check_and_repair_state(BE)

        !--------------------------------------------------------------------
        ! Prepare all fields FABM needs to compute source terms (e.g., light)
        !---------------------------------------------------------------------
        call BE%model%prepare_inputs(istep_rk) !Providing the main time-step number as set_domain received dt_main (t, year, month, day, seconds)

        !-------------------------------------------------------------
        ! Obtaining tendencies (dC/dt) and surface and bottom fluxes
        !-------------------------------------------------------------
        BE%flux_sf = 0.0_rk; BE%tendency_sf  = 0.0_rk                 ! Zeroing before receiving surface fluxes and sources
        call BE%model%get_surface_sources(BE%flux_sf, BE%tendency_sf) ! Retrieve surface fluxes and sources from FABM

        BE%flux_bt = 0.0_rk; BE%tendency_bt  = 0.0_rk                 ! Zero before receiving bottom fluxes and sources
        call BE%model%get_bottom_sources(BE%flux_bt, BE%tendency_bt)  ! Retrieve bottom fluxes and sources from FABM

        BE%tendency_int = 0.0_rk
        call BE%model%get_interior_sources(1, nz, BE%tendency_int)

        ! Include contribution of bottom and surface fluxes to tendencies at bottom and surface
        if (nint > 0) then
            do ivar = 1, nint
                ! Bottom flux: positive is flux into the column, negative out of the column
                ! Flux = mass/(Area*dt)
                ! Concentration change due to flux = mass flux/volume = (Flux Area)/(Area dz) = Flux/dz
                BE%tendency_int(1,ivar)  = BE%tendency_int(1,ivar)  + BE%flux_bt(ivar) / BE%grid%dz(1)

                ! Surface flux: positive is flux into the water, negative out of the water
                ! Concentration change due to flux = mass flux/volume = (Flux Area)/(Area dz) = Flux/dz
                BE%tendency_int(nz,ivar) = BE%tendency_int(nz,ivar) + BE%flux_sf(ivar) / BE%grid%dz(nz)
            end do
        end if    

        !------------------------------------------------
        ! Let FABM compute any remaining outputs (diagnostics) 
        !------------------------------------------------
        call BE%model%finalize_outputs()
       
        !------------------------------------------------------------
        ! Apply vertical residual movement due to sinking or floating 
        !------------------------------------------------------------
        if (nint>0) then
            BE%velocity = 0.0_rk
            call BE%model%get_vertical_movement(1,nz, BE%velocity)   ! Obtain vertical velocities from FABM
            do ivar=1, nint
                call apply_vertical_transport(BE%BS%interior_state(:,ivar), BE%grid,&
                                              vel_center=BE%velocity(:,ivar), dt=dt_main)
            end do
        end if
        
        !-------------------------------------------------------------------------
        ! Mix internal tracers vertically due to turbulent diffusion.
        !--------------------------------------------------------------------------
        if (nint>0) then
            do ivar=1, nint
                ! Only internal mixing is applied: no diffusive flux enters or leaves through surface or bottom
                call scalar_diffusion(Var=BE%BS%interior_state(:,ivar), N=nz, dt=dt_main, h=BE%grid%dz, &
                                      Kz=BE%BS%vert_diff, cnpar=BE%params%cnpar, tricoef=BE%trid, ierr=ierr)
            end do
        end if

        !-----------------------------------------------------------------------------
        ! Repir state for interior tracers after vertical redistribution of tracers
        !----------------------------------------------------------------------------
        call check_and_repair_state(BE)
    

        ! ----------------------------------------------------------------------------
        ! Compute number of substeps to guarantee numerical stability for intergation of tendencies
        !-----------------------------------------------------------------------------
        call compute_bio_substeps(BE, BE%tendency_int, BE%tendency_sf, BE%tendency_bt, dt_main, BE%params%frac_max, nsub, dt_sub)

        !-----------------------------------------------------------------------------
        ! Integrate tracers using exflicit Forward Euler
        !-----------------------------------------------------------------------------    
        do isub = 1, nsub
        ! Subcycle to achieve stability (Euler can be stiff)
            ! Interior
            if (nint>0) then
                do ivar = 1, nint
                    do k = 1, nz
                        BE%BS%interior_state(k,ivar) = BE%BS%interior_state(k,ivar) + dt_sub * BE%tendency_int(k,ivar)
                    end do
                end do
            end if

            ! Surface
            if (nsfc>0) then
                do ivar = 1, nsfc
                    BE%BS%surface_state(ivar) = BE%BS%surface_state(ivar) + dt_sub * BE%tendency_sf(ivar)
                end do
            end if 

            ! Bottom
            if (nbtm>0) then
                do ivar = 1, nbtm
                    BE%BS%bottom_state(ivar)  = BE%BS%bottom_state(ivar)  + dt_sub * BE%tendency_bt(ivar)
                end do
            end if
        end do      

        ! Repair state, if needed, to let FABM restore all state variables within their registered bounds.
        call check_and_repair_state(BE)      

!!!
    write(*,*) 'nsub=', nsub, ' dt_sub=', dt_sub, ' valid_int=', BE%valid_int
    write(*,'(/,A,I0)') '--- Bio profile at main step = ', istep_main
    write(*,'(A)') '    Depth(m)   Tracers (one column per variable)'

    ! Optional: print tracer names (truncated to 10 chars)
    write(*,'(A)', advance='no') '    z[m]     '
    do ivar = 1, nint
        write(*,'(1X,A10)', advance='no') trim(BE%BS%intvar_names(ivar)(1:min(10,len_trim(BE%BS%intvar_names(ivar)))))
    end do
    write(*,*)
    write(*,'(F9.3,1X,1P,100E12.3)') BE%grid%z(nz), (BE%flux_sf(ivar), ivar=1,nint)
    write(*,'(F9.3,1X,1P,100E12.3)') BE%grid%z(1), (BE%flux_bt(ivar), ivar=1,nint)
    write(*,'(F9.3,1X,1P,100E12.3)') BE%grid%z(1), (BE%velocity(1,ivar), ivar=1,nint)

    ! Values: one row per depth
    !do k = 1, nz
    !    write(*,'(F9.3,1X,1P,100E12.3)') BE%grid%z(k), (BE%BS%interior_state(k,ivar), ivar=1,nint)
    !end do
!!!
    end subroutine integrate_bio_fabm

    !==============================================================
    ! Clears memory associated with this module
    !==============================================================
    subroutine end_bio_fabm(BE)
        type(BioEnv), intent(inout) :: BE


        ! Report number of times that variables where repaired
        if (BE%params%repair) then
            if (BE%BS%n_interior > 0) write(*,'(A,I0,A)') 'FABM repaired the interior variables ', BE%nrepair_int, ' time(s).'
            if (BE%BS%n_surface > 0)  write(*,'(A,I0,A)') 'FABM repaired the surface variables ', BE%nrepair_sfc, ' time(s).'
            if (BE%BS%n_bottom > 0)   write(*,'(A,I0,A)') 'FABM repaired the bottom variables ', BE%nrepair_btm, ' time(s).'
        end if

        !-------------------------------------------------
        ! Clear FABM model instance 
        !-------------------------------------------------
        if (associated(BE%model)) then
            call BE%model%finalize()
            nullify(BE%model)
        end if

        !-------------------------------------------------
        ! Deallocate BioState arrays
        !-------------------------------------------------
        if (allocated(BE%BS%interior_state)) deallocate(BE%BS%interior_state)
        if (allocated(BE%BS%bottom_state))   deallocate(BE%BS%bottom_state)
        if (allocated(BE%BS%surface_state))  deallocate(BE%BS%surface_state)

        if (allocated(BE%BS%temp))      deallocate(BE%BS%temp)
        if (allocated(BE%BS%sal))       deallocate(BE%BS%sal)
        if (allocated(BE%BS%rho))       deallocate(BE%BS%rho)
        if (allocated(BE%BS%pres))      deallocate(BE%BS%pres)
        if (allocated(BE%BS%swr))       deallocate(BE%BS%swr)
        if (allocated(BE%BS%par))       deallocate(BE%BS%par)
        if (allocated(BE%BS%vert_diff)) deallocate(BE%BS%vert_diff)

        ! Variable-name arrays
        if (allocated(BE%BS%intvar_names)) deallocate(BE%BS%intvar_names)
        if (allocated(BE%BS%sfcvar_names)) deallocate(BE%BS%sfcvar_names)
        if (allocated(BE%BS%btmvar_names)) deallocate(BE%BS%btmvar_names)

        !-------------------------------------------------
        ! Deallocate BioEnv working arrays
        !-------------------------------------------------
        if (allocated(BE%velocity))     deallocate(BE%velocity)
        if (allocated(BE%tendency_int)) deallocate(BE%tendency_int)
        if (allocated(BE%tendency_sf))  deallocate(BE%tendency_sf)
        if (allocated(BE%tendency_bt))  deallocate(BE%tendency_bt)
        if (allocated(BE%flux_sf))      deallocate(BE%flux_sf)
        if (allocated(BE%flux_bt))      deallocate(BE%flux_bt)

        ! Clear Tridiagonal workspace 
        call clear_tridiag(BE%trid)

        !-------------------------------------------------
        ! Reset counters and flags
        !-------------------------------------------------
        BE%BS%n_interior = 0
        BE%BS%n_surface  = 0
        BE%BS%n_bottom   = 0
        BE%BS%n_total    = 0

        BE%need_temp    = .false.
        BE%need_salt    = .false.
        BE%need_rho     = .false.
        BE%need_pres    = .false.
        BE%need_par     = .false.
        BE%need_swr     = .false.
        BE%need_windspd = .false.
        BE%need_slp     = .false.
        BE%need_par_sfc = .false.
        BE%need_swr_sfc = .false.
        BE%need_cloud   = .false.
        BE%need_stressb = .false.
        BE%need_co2     = .false.

        BE%is_init = .false.
    end subroutine end_bio_fabm



    !============== INTERNAL HELPERS ===============================

    ! Provides environment data to FABM
    subroutine link_environment_data(PS, FS, BE)
        type(PhysicsState),       intent(in) :: PS
        type(ForcingSnapshot),    intent(in) :: FS
        type(BioEnv),          intent(inout) :: BE
 !! When sediments plugged, allocate first with the size of the full grid       
        BE%BS%temp = PS%temp  !Later decide what to do with the sediment layers
        BE%BS%sal  = PS%sal
        BE%BS%rho  = PS%rho
        Be%BS%vert_diff = PS%Kz + mol_diff  ! Adding molecular diffusivity
        BE%BS%short_rad = FS%short_rad
        BE%BS%wind_spd  = FS%wind_spd
        BE%BS%slp       = FS%slp * 100.0_rk ! Converting from hPa to Pa
        BE%BS%co2_air   = FS%co2_air
        BE%BS%stressb   = PS%stressb
        BE%BS%par_sfc   = BE%BS%short_rad * 0.45_rk ! Multiply instead by par_fraction as a parameter

        ! Passing all these variables even if they are not needed as the information is already available
        call BE%model%link_interior_data(BE%id_temp, BE%BS%temp)
        call BE%model%link_interior_data(BE%id_salt, BE%BS%sal)
        call BE%model%link_interior_data(BE%id_rho, BE%BS%rho)
        if (BE%need_pres) then
            call compute_pressure(BE%grid%dz, BE%BS%rho, BE%BS%pres)
            call BE%model%link_interior_data(BE%id_pres, BE%BS%pres)
        end if
        ! Profile of shortwave radiation
        if(BE%need_swr) then
            call compute_SW_profile(BE%wat_grid%dz, BE%BS%short_rad, BE%BS%swr)  ! CORRECT when sediments added
            call BE%model%link_interior_data(BE%id_swr, BE%BS%swr)
        end if
        ! PAR Profile
        if(BE%need_par) then
            call compute_PAR_profile(BE%wat_grid%dz, BE%BS%short_rad, BE%BS%par)   ! CORRECT when sediments are added
            call BE%model%link_interior_data(BE%id_par, BE%BS%par)
        end if
        ! Linking surface data
        call  BE%model%link_horizontal_data(BE%id_windspd, BE%BS%wind_spd)
        call  BE%model%link_horizontal_data(BE%id_slp, BE%BS%slp)
        call  BE%model%link_horizontal_data(BE%id_swr_sfc, BE%BS%short_rad)
        call  BE%model%link_horizontal_data(BE%id_co2, BE%BS%co2_air)
        call  BE%model%link_horizontal_data(BE%id_par_sfc, BE%BS%par_sfc)
        if(BE%need_stressb) call  BE%model%link_horizontal_data(BE%id_stressb, BE%BS%stressb)

    end subroutine link_environment_data 

    subroutine check_and_repair_state(BE)
        type(BioEnv), intent(inout) :: BE

        logical :: repair
        repair = BE%params%repair

        ! Repair state, if needed, to let FABM restore all state variables within their registered bounds.
        if (BE%BS%n_interior > 0) then
             call BE%model%check_interior_state(1, BE%grid%nz, repair, BE%valid_int)
             if (.not. BE%valid_int .and. repair) BE%nrepair_int = BE%nrepair_int + 1 ! Update counters
        end if
        if (BE%BS%n_surface > 0) then
            call BE%model%check_surface_state(repair, BE%valid_sfc)
             if (.not. BE%valid_sfc .and. repair) BE%nrepair_sfc = BE%nrepair_sfc + 1 ! Update counters
        end if
        if (BE%BS%n_bottom > 0) then
            call BE%model%check_bottom_state(repair, BE%valid_btm)
            if (.not. BE%valid_btm .and. repair) BE%nrepair_btm = BE%nrepair_btm + 1 ! Update counters
        end if 

        if (.not. (BE%valid_int .and. BE%valid_sfc .and. BE%valid_btm)) then
            if (.not. repair) stop 'integrate_bio_fabm: invalid state found in a biogeochemistry variable. '
        end if    
    end subroutine check_and_repair_state

end module bio_main