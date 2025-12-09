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
    use time_types,          only: DateTime
    use tridiagonal,         only: init_tridiag, clear_tridiag
    use vertical_transport,  only: apply_vertical_transport, velocity_at_interfaces
    use vertical_mixing,     only: scalar_diffusion
    use variable_registry,   only: register_variable, output_all_variables


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

        integer  :: ivar, nz, nint, nsfc, nbtm
        integer :: n_total_int_diag, n_save_int, j
        integer :: n_total_hz_diag, n_save_hz
        real(rk) :: dt_main

        write(*,'(A)') 'Initialising biogeochemistry via FABM...'

        BE%grid     = grid   ! Full grid
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
            do ivar = 1, nint
                ! Point FABM to the current concentration state for biogeochemical tracers
                call BE%model%link_interior_state_data(ivar, BE%BS%interior_state(:,ivar))
                ! Register metadata for variables
                call register_variable(BE%int_vars, name=trim(BE%model%interior_state_variables(ivar)%name),  &
                                       long_name=trim(BE%model%interior_state_variables(ivar)%long_name),        &
                                       units=trim(BE%model%interior_state_variables(ivar)%units),               &
                                       minimum=BE%model%interior_state_variables(ivar)%minimum,                 &
                                       maximum=BE%model%interior_state_variables(ivar)%maximum,                 &
                                       missing_value=BE%model%interior_state_variables(ivar)%missing_value,     &
                                       vert_coord='centre', n_space_dims=1, data_1d=BE%BS%interior_state(:,ivar))   
            end do
            if (allocated(BE%int_vars)) call output_all_variables(BE%int_vars)
        end if

        ! Allocating working arrays        
        allocate(BE%tendency_int(nz,nint))                        ! Sources (dC/dt) 
        allocate(BE%velocity(nz,nint))                            ! Velocities
        allocate(BE%vel_faces(0:nz,nint))                         ! Velocities at interfaces
        allocate(BE%flux_sf(nint));    allocate(BE%flux_bt(nint)) ! Fluxes at the surface and bottom of interior variables 
        allocate(BE%BS%vert_diff(0:nz))                           ! Vertical diffusivities

        ! --- Bottom-only state variables ---
        nbtm = size(BE%model%bottom_state_variables)
        BE%BS%n_bottom = nbtm
        if (nbtm > 0) then
            if (allocated(BE%BS%bottom_state)) deallocate(BE%BS%bottom_state)
            allocate(BE%BS%bottom_state(nbtm))            
            ! Link FABM to the bottom state vector
            call BE%model%link_all_bottom_state_data(BE%BS%bottom_state)
            do ivar=1, nbtm
                call register_variable(BE%btm_vars, name=trim(BE%model%bottom_state_variables(ivar)%name),  &
                                       long_name=trim(BE%model%bottom_state_variables(ivar)%long_name),        &
                                       units=trim(BE%model%bottom_state_variables(ivar)%units),               &
                                       minimum=BE%model%bottom_state_variables(ivar)%minimum,                 &
                                       maximum=BE%model%bottom_state_variables(ivar)%maximum,                 &
                                       missing_value=BE%model%bottom_state_variables(ivar)%missing_value,     &
                                       vert_coord='bottom', n_space_dims=0, data_0d=BE%BS%bottom_state(ivar)) 
            end do
            if (allocated(BE%btm_vars)) call output_all_variables(BE%btm_vars)
        end if
        ! Allocate working array for sources (dC/dt) at the bottom
        allocate(BE%tendency_bt(nbtm))

        ! --- Surface-only state variables ---
        nsfc = size(BE%model%surface_state_variables)
        BE%BS%n_surface = nsfc
        if (nsfc > 0) then
            if (allocated(BE%BS%surface_state)) deallocate(BE%BS%surface_state)
            allocate(BE%BS%surface_state(nsfc))         
            ! Link FABM to the surface state vector
            call BE%model%link_all_surface_state_data(BE%BS%surface_state)
            do ivar=1, nsfc
                call register_variable(BE%sfc_vars, name=trim(BE%model%surface_state_variables(ivar)%name),  &
                                       long_name=trim(BE%model%surface_state_variables(ivar)%long_name),        &
                                       units=trim(BE%model%surface_state_variables(ivar)%units),               &
                                       minimum=BE%model%surface_state_variables(ivar)%minimum,                 &
                                       maximum=BE%model%surface_state_variables(ivar)%maximum,                 &
                                       missing_value=BE%model%surface_state_variables(ivar)%missing_value,     &
                                       vert_coord='surface', n_space_dims=0, data_0d=BE%BS%surface_state(ivar)) 
            end do
            if (allocated(BE%sfc_vars)) call output_all_variables(BE%sfc_vars)
        end if
        ! Allocate working array for sources (dC/dt) at the surface
        allocate(BE%tendency_sf(nsfc));                

        
        BE%BS%n_total = nint + nsfc + nbtm  ! Total number of state variables

        !---- Interior diagnostic variables -----------
        n_total_int_diag = size(BE%model%interior_diagnostic_variables)
        n_save_int       = 0
        ! Count diagnostics that bio models marked as save=true
        do ivar = 1, n_total_int_diag
            if (BE%model%interior_diagnostic_variables(ivar)%save) n_save_int = n_save_int + 1
        end do
        BE%n_diag_int = n_save_int

        if (n_save_int > 0) then
            if (allocated(BE%diag_int))       deallocate(BE%diag_int)
            if (allocated(BE%diag_int_index)) deallocate(BE%diag_int_index)
            if (allocated(BE%diag_int_vars))  deallocate(BE%diag_int_vars)

            allocate(BE%diag_int(nz, n_save_int))
            allocate(BE%diag_int_index(n_save_int))

            j = 0
            do ivar = 1, n_total_int_diag
                if (.not. BE%model%interior_diagnostic_variables(ivar)%save) cycle
                j = j + 1
                BE%diag_int_index(j) = ivar   ! map our j-th diagnostic to FABM index i

                call register_variable( BE%diag_int_vars,                                    &
                                        name        = trim(BE%model%interior_diagnostic_variables(ivar)%name),      &
                                        long_name   = trim(BE%model%interior_diagnostic_variables(ivar)%long_name), &
                                        units       = trim(BE%model%interior_diagnostic_variables(ivar)%units),     &
                                        minimum=BE%model%interior_diagnostic_variables(ivar)%minimum,                 &
                                        maximum=BE%model%interior_diagnostic_variables(ivar)%maximum,                 &
                                        missing_value=BE%model%interior_diagnostic_variables(ivar)%missing_value,     &
                                        vert_coord  = 'centre',                                                 &
                                        n_space_dims= 1,                                                        &
                                        data_1d     = BE%diag_int(:, j) )
            end do
            if (allocated(BE%diag_int_vars)) call output_all_variables(BE%diag_int_vars)
        end if

        ! Horizontal diagnostics with save=true
        n_total_hz_diag = size(BE%model%horizontal_diagnostic_variables)
        n_save_hz       = 0

        do ivar = 1, n_total_hz_diag
            if (BE%model%horizontal_diagnostic_variables(ivar)%save) n_save_hz = n_save_hz + 1
        end do
        BE%n_diag_hz = n_save_hz

        if (n_save_hz > 0) then
            if (allocated(BE%diag_hz))       deallocate(BE%diag_hz)
            if (allocated(BE%diag_hz_index)) deallocate(BE%diag_hz_index)
            if (allocated(BE%diag_hz_vars))  deallocate(BE%diag_hz_vars)

            allocate(BE%diag_hz(n_save_hz))
            allocate(BE%diag_hz_index(n_save_hz))

            j = 0
            do ivar = 1, n_total_hz_diag
                if (.not. BE%model%horizontal_diagnostic_variables(ivar)%save) cycle
                j = j + 1
                BE%diag_hz_index(j) = ivar   ! map our j-th diagnostic to FABM index

                call register_variable( BE%diag_hz_vars,                                         &
                                        name        = trim(BE%model%horizontal_diagnostic_variables(ivar)%name),        &
                                        long_name   = trim(BE%model%horizontal_diagnostic_variables(ivar)%long_name),   &
                                        units       = trim(BE%model%horizontal_diagnostic_variables(ivar)%units),       &
                                        minimum=BE%model%horizontal_diagnostic_variables(ivar)%minimum,                 &
                                        maximum=BE%model%horizontal_diagnostic_variables(ivar)%maximum,                 &
                                        missing_value=BE%model%horizontal_diagnostic_variables(ivar)%missing_value,     &
                                        vert_coord  = 'none',                                                        &
                                        n_space_dims= 0,                                                             &
                                        data_0d     = BE%diag_hz(j))
            end do
            if (allocated(BE%diag_hz_vars)) call output_all_variables(BE%diag_hz_vars)
        end if

        ! For conserved Quantities
        if (BE%params%output_conserved) then
            BE%n_conserved = size(BE%model%conserved_quantities)
            if (BE%n_conserved > 0) then
                allocate(BE%conserved_interior(nz, BE%n_conserved))
                allocate(BE%conserved_boundary(BE%n_conserved))
                allocate(BE%conserved_total(BE%n_conserved))

                BE%conserved_interior = 0.0_rk
                BE%conserved_boundary = 0.0_rk
                BE%conserved_total    = 0.0_rk
                ! register variables for output      
                do ivar=1, BE%n_conserved
                    call register_variable(BE%conserved_vars,                                 &
                                            name='conserved_total_'//trim(BE%model%conserved_quantities(ivar)%name),  &
                                            long_name='Column-integrated total of '//trim(BE%model%conserved_quantities(ivar)%long_name),        &
                                            units=trim(BE%model%conserved_quantities(ivar)%units)//'*m',         &
                                            minimum=BE%model%conserved_quantities(ivar)%minimum,                 &
                                            maximum=BE%model%conserved_quantities(ivar)%maximum,                 &
                                            missing_value=BE%model%conserved_quantities(ivar)%missing_value,     &
                                            vert_coord='surface', n_space_dims=0, data_0d=BE%conserved_total(ivar)) 
                end do
                if (allocated(BE%conserved_vars)) call output_all_variables(BE%conserved_vars)
            end if
        end if

        ! Provide location info to FABM
        call BE%model%link_interior_data(fabm_standard_variables%depth,BE%grid%z)                ! layer depths at centres
        call BE%model%link_interior_data(fabm_standard_variables%cell_thickness,BE%grid%dz)      ! layer thicknesses
        call BE%model%link_horizontal_data(fabm_standard_variables%latitude, location%lat)       ! Latitude
        call BE%model%link_horizontal_data(fabm_standard_variables%longitude, location%lon)      ! Longitude
        call BE%model%link_horizontal_data(fabm_standard_variables%bottom_depth, location%depth) ! Bottom depth
        ! call model%set_bottom_index(bottom_indices)

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
        ! Scalar
        BE%id_yearday = BE%model%get_scalar_variable_id(fabm_standard_variables%number_of_days_since_start_of_the_year)
        call BE%model%link_scalar(BE%id_yearday, BE%BS%doy)
        
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

    !!! Change when sediments are integrated   
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

        ! Initialising counters
        BE%nrepair_int = 0; BE%nrepair_btm=0; BE%nrepair_sfc=0;
        BE%valid_int = .true.; BE%valid_btm=.true.; BE%valid_sfc=.true.
        
        BE%is_init = .true.
        write(*,'("✓ Biogeochemistry initialised successfully with ",I0," state variables:")') BE%BS%n_total
        write(*,'(" ",I0," interior, ",I0," surface and ",I0," bottom variables")') BE%BS%n_interior, BE%BS%n_surface, BE%BS%n_bottom
        write(*,'(" ",I0," interior and ",I0," horizontal diagnostic variable(s).")') &
        BE%n_diag_int, BE%n_diag_hz
        if (BE%params%output_conserved) then
            if (BE%n_conserved > 0) then
                write(*,'(" ",I0," conserved quantity(ies) will be tracked as column-integrated totals.")') &
                    BE%n_conserved
            else
                write(*,'(" FABM reports no conserved quantities for this configuration.")')
            end if
        else
            if (size(BE%model%conserved_quantities) > 0) then
                write(*,'(" Conserved quantities are available (",I0,") but output is disabled (output_conserved_qt: no).")') &
                    size(BE%model%conserved_quantities)
            end if
        end if
    end subroutine init_bio_fabm

    !=====================================================================================
    ! Update biogeochemistry for one timestep.
    !   - Communicates to FABM current environmental state.
    !   - Gets from FABM source terms and boundary fluxes.
    !   - Applies vertical mixing and moves tracers according to reported velocities.
    !   - Integrates tracer tendencies (sources: dC/dt) subcycling if needed for stability.
    !=====================================================================================
    subroutine integrate_bio_fabm(BE, PS, FS, timestep, istep_main, date, sec_of_day, doy_fraction)
        type(BioEnv),          intent(inout) :: BE
        type(PhysicsState),    intent(in)    :: PS
        type(ForcingSnapshot), intent(in)    :: FS
        integer(lk),           intent(in)    :: timestep
        integer(lk),           intent(in)    :: istep_main
        type(DateTime),        intent(in)    :: date
        real(rk),              intent(in)    :: sec_of_day
        real(rk),              intent(in)    :: doy_fraction   ! 0-based day of year + fraction


        integer  :: nz, ivar, k, nint, nsfc, nbtm
        real(rk) :: istep_rk, dt_main, dt_sub
        integer  :: nsub, isub, nsub_adv, nsub_diff, nsub_bio
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

        BE%BS%doy = doy_fraction

        !---------------------------------------------
        ! Update FABM with environment data
        !--------------------------------------------
        call link_environment_data(PS, FS, BE)

        ! Repair state, if needed, to let FABM restore all state variables within their registered bounds.
        call check_and_repair_state(BE)

        !--------------------------------------------------------------------
        ! Prepare all fields FABM needs to compute source terms (e.g., light)
        !---------------------------------------------------------------------
        call BE%model%prepare_inputs(istep_rk, date%year, date%month, date%day, sec_of_day) !Providing the main time-step number as set_domain received dt_main 

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

        call update_bio_diagnostics(BE)        

        BE%velocity = 0.0_rk
        BE%vel_faces = 0._rk   
        call BE%model%get_vertical_movement(1,nz, BE%velocity)   ! Obtain vertical velocities from FABM
        call velocity_at_interfaces(BE%velocity, BE%grid, BE%vel_faces)  ! Caluclate velocities at layer interfaces

        ! Compute number of subcycles needed for numerical stabiluty
        call compute_bio_substeps(BE%vel_faces, BE%grid%dz, Kz=BE%BS%vert_diff,         &
                                  cnpar=BE%params%cnpar, BE=BE, dt_main=dt_main,        &
                                  frac_max=BE%params%frac_max, dt_min=BE%params%min_dt, &
                                  nsub_adv=nsub_adv, nsub_diff=nsub_diff,               &
                                  nsub_bio=nsub_bio, nsub=nsub, dt_sub=dt_sub)

        ! Inner loop to maintain numerical stabilty
        do isub=1, nsub
            !------------------------------------------------------------
            ! Apply vertical residual movement and mix tracers vertically
            !------------------------------------------------------------
            if (nint>0) then
                
                do ivar=1, nint
                    if (any(BE%vel_faces(:, ivar) /= 0._rk)) then  
                        ! Only if the tracer has vertical motion
                        ! Move tracers due to sinking or floating 
                        call apply_vertical_transport(BE%BS%interior_state(:,ivar), BE%grid,&
                                                    w_face=BE%vel_faces(:,ivar), dt=dt_sub)
                    end if
                    ! Mix internal tracers vertically due to turbulent diffusion.
                    call scalar_diffusion(Var=BE%BS%interior_state(:,ivar), N=nz, dt=dt_sub, h=BE%grid%dz, &
                                          Kz=BE%BS%vert_diff, cnpar=BE%params%cnpar, tricoef=BE%trid, ierr=ierr)
                end do
            end if   
            
            call BE%model%link_all_interior_state_data(BE%BS%interior_state)

            !-----------------------------------------------------------------------------
            ! Repir state for interior tracers after vertical redistribution of tracers
            !----------------------------------------------------------------------------
            call check_and_repair_state(BE)
            if (.not. BE%valid_int) then
                write(*,*) "After transport"
                stop 1
            end if
            
            !-----------------------------------------------------------------------------
            ! Integrate tracers using exflicit Forward Euler
            !-----------------------------------------------------------------------------    
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
            ! Repair state, if needed, to let FABM restore all state variables within their registered bounds.
            call check_and_repair_state(BE)       
        end do   
        
        ! Compute totals for conserved quantities if needed
        if (BE%params%output_conserved .and. BE%n_conserved > 0) then
            call BE%model%get_interior_conserved_quantities(1, nz, BE%conserved_interior)
            call BE%model%get_horizontal_conserved_quantities(BE%conserved_boundary)
            do ivar = 1, BE%n_conserved
                BE%conserved_total(ivar) = sum(BE%grid%dz(1:nz) * BE%conserved_interior(1:nz, ivar)) + BE%conserved_boundary(ivar)
            end do
        end if           

    end subroutine integrate_bio_fabm

    !==============================================================
    ! Clears memory associated with bio_main
    !==============================================================
    subroutine end_bio_fabm(BE)
        type(BioEnv), intent(inout) :: BE

        integer :: i

        ! Report number of times that variables where repaired
        if (BE%params%repair) then
            if (BE%nrepair_int>0 .or. BE%nrepair_sfc>0 .or. BE%nrepair_btm>0) then
                write(*,*) 'Warning: FABM repaired some state variables;'
                write(*,*) '         reducing the time step may prevent this.'
            end if
            if (BE%BS%n_interior > 0) write(*,'(A,I0,A)') 'FABM repaired the interior variables ', BE%nrepair_int, ' time(s).'
            if (BE%BS%n_surface > 0)  write(*,'(A,I0,A)') 'FABM repaired the surface variables ', BE%nrepair_sfc, ' time(s).'
            if (BE%BS%n_bottom > 0)   write(*,'(A,I0,A)') 'FABM repaired the bottom variables ', BE%nrepair_btm, ' time(s).'
        end if

        ! Clear FABM model instance 
        if (associated(BE%model)) then
            call BE%model%finalize()
            nullify(BE%model)
        end if

        ! Deallocate BioState arrays
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

        !Clean variable Metadata
        ! Interior bio variables
        if (allocated(BE%int_vars)) then
            do i = 1, size(BE%int_vars)
                nullify(BE%int_vars(i)%data_0d)
                nullify(BE%int_vars(i)%data_1d)
            end do
            deallocate(BE%int_vars)
        end if

        ! Surface bio variables
        if (allocated(BE%sfc_vars)) then
            do i = 1, size(BE%sfc_vars)
                nullify(BE%sfc_vars(i)%data_0d)
                nullify(BE%sfc_vars(i)%data_1d)
            end do
            deallocate(BE%sfc_vars)
        end if

        ! Bottom bio variables
        if (allocated(BE%btm_vars)) then
            do i = 1, size(BE%btm_vars)
                nullify(BE%btm_vars(i)%data_0d)
                nullify(BE%btm_vars(i)%data_1d)
            end do
            deallocate(BE%btm_vars)
        end if      
        
        ! Diagnostics interior and horizontal
        if (allocated(BE%diag_int))       deallocate(BE%diag_int)
        if (allocated(BE%diag_int_index)) deallocate(BE%diag_int_index)
        if (allocated(BE%diag_hz))        deallocate(BE%diag_hz)
        if (allocated(BE%diag_hz_index))  deallocate(BE%diag_hz_index)
        ! Clearing metadata for interior diagnostics
        if (allocated(BE%diag_int_vars)) then
            do i = 1, size(BE%diag_int_vars)
                nullify(BE%diag_int_vars(i)%data_0d)
                nullify(BE%diag_int_vars(i)%data_1d)
            end do
            deallocate(BE%diag_int_vars)
        end if
        ! Clearing metadata for horizontale diagnostics
        if (allocated(BE%diag_hz_vars)) then
            do i = 1, size(BE%diag_hz_vars)
                nullify(BE%diag_hz_vars(i)%data_0d)
                nullify(BE%diag_hz_vars(i)%data_1d)
            end do
            deallocate(BE%diag_hz_vars)
        end if

        BE%n_diag_int = 0
        BE%n_diag_hz  = 0

        ! Deallocate arrays for conserved quantities
        if (allocated(BE%conserved_interior)) deallocate(BE%conserved_interior)
        if (allocated(BE%conserved_boundary)) deallocate(BE%conserved_boundary)
        if (allocated(BE%conserved_total))    deallocate(BE%conserved_total)

        if (allocated(BE%conserved_vars)) then
            do i = 1, size(BE%conserved_vars)
                nullify(BE%conserved_vars(i)%data_0d)
                nullify(BE%conserved_vars(i)%data_1d)
            end do
            deallocate(BE%conserved_vars)
        end if
        BE%n_conserved = 0


        ! Deallocate BioEnv working arrays
        if (allocated(BE%velocity))     deallocate(BE%velocity)
        if (allocated(BE%tendency_int)) deallocate(BE%tendency_int)
        if (allocated(BE%tendency_sf))  deallocate(BE%tendency_sf)
        if (allocated(BE%tendency_bt))  deallocate(BE%tendency_bt)
        if (allocated(BE%flux_sf))      deallocate(BE%flux_sf)
        if (allocated(BE%flux_bt))      deallocate(BE%flux_bt)

        ! Clear Tridiagonal workspace 
        call clear_tridiag(BE%trid)

        ! Reset counters and flags
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
        BE%BS%vert_diff = PS%Kz * 0.05_rk !+ mol_diff  ! Adding molecular diffusivity
        BE%BS%short_rad = FS%short_rad
        BE%BS%wind_spd  = PS%wind_speed
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
            call compute_SW_profile(BE%grid%dz, BE%BS%short_rad, BE%BS%swr)  ! CORRECT when sediments added
            call BE%model%link_interior_data(BE%id_swr, BE%BS%swr)
        end if
        ! PAR Profile
        if(BE%need_par) then
            call compute_PAR_profile(BE%grid%dz, BE%BS%short_rad, BE%BS%par)   ! CORRECT when sediments are added
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

        if (.not. (BE%valid_int .and. BE%valid_sfc .and. BE%valid_btm) .and. .not. repair) then            
            write(*,*) 'ERROR in integrate_bio_fabm: invalid biogeochemical state detected.'
            write(*,*) '      Try reducing the main time-step and/or the minimum substep (min_timestep).'
            stop 1
        end if  
    end subroutine check_and_repair_state

    subroutine update_bio_diagnostics(BE)
        type(BioEnv), intent(inout) :: BE

        integer :: j, idx

        ! Interior diagnostics 
        if (BE%n_diag_int > 0) then
            do j = 1, BE%n_diag_int
                idx   = BE%diag_int_index(j)    ! FABM index
                BE%diag_int(:, j) = BE%model%get_interior_diagnostic_data(idx)
            end do
        end if

        ! Horizontal diagnostics 
        if (BE%n_diag_hz > 0) then
            do j = 1, BE%n_diag_hz
                idx   = BE%diag_hz_index(j)
                BE%diag_hz(j) = BE%model%get_horizontal_diagnostic_data(idx)
            end do
        end if
    end subroutine update_bio_diagnostics
end module bio_main