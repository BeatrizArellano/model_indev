module bio_main
    use fabm
    use bio_types,        only: BioEnv
    use geo_utils,        only: LocationInfo
    use physics_types,    only: PhysicsState
    use forcing_manager,  only: ForcingSnapshot
    use grids,            only: VerticalGrid   
    use light,            only: compute_SW_profile, compute_PAR_profile
    use pressure,         only: compute_pressure 
    use precision_types,  only: rk, lk


  implicit none
  private

  public :: BioEnv, init_bio_fabm


contains

    subroutine init_bio_fabm(location, grid, timestep, PS, FS, BE)
        type(LocationInfo),       intent(in) :: location
        type(VerticalGrid),       intent(in) :: grid
        integer(lk), intent(in)              :: timestep
        type(PhysicsState),       intent(in) :: PS
        type(ForcingSnapshot),    intent(in) :: FS
        type(BioEnv),          intent(inout) :: BE

        integer  :: maxlen, ivar
        real(rk) :: dt_main

        write(*,'(A)') 'Initialising biogeochemistry via FABM...'

        BE%grid     = grid   ! Full grid
        BE%wat_grid = grid   ! For now they're the same grid
        dt_main     = real(timestep, rk)

        ! Initialize (reads FABM configuration from fabm.yaml)
        ! After this the number of biogeochemical variables is fixed.
        BE%model => fabm_create_model()
    

        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        call BE%model%set_domain(BE%grid%nz, seconds_per_time_unit=dt_main)

        ! ---- Interior state variables ---------------
        BE%BS%n_interior = size(BE%model%interior_state_variables)
        if (BE%BS%n_interior > 0) then
            ! Allocate the array for interior tracers concentrations
            if(allocated(BE%BS%interior_state)) deallocate(BE%BS%interior_state)
            allocate(BE%BS%interior_state(BE%grid%nz, BE%BS%n_interior))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%interior_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%intvar_names(BE%BS%n_interior))
            ! Point FABM to the current concentration state for biogeochemical tracers
            do ivar = 1, BE%BS%n_interior
                BE%BS%intvar_names(ivar) = trim(BE%model%interior_state_variables(ivar)%name) ! Retrieve names of tracers
                call BE%model%link_interior_state_data(ivar, BE%BS%interior_state(:,ivar))
            end do
        end if


        ! --- Bottom-only state variables ---
        BE%BS%n_bottom = size(BE%model%bottom_state_variables)
        if (BE%BS%n_bottom > 0) then
            if (allocated(BE%BS%bottom_state)) deallocate(BE%BS%bottom_state)
            allocate(BE%BS%bottom_state(BE%BS%n_bottom))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%bottom_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%btmvar_names(BE%BS%n_bottom))
            ! Link FABM to the bottom state vector
            call BE%model%link_all_bottom_state_data(BE%BS%bottom_state)
            do ivar=1, BE%BS%n_bottom
                BE%BS%btmvar_names(ivar) = trim(BE%model%bottom_state_variables(ivar)%name) ! Retrieve names of tracers
            end do
        end if

        ! --- Surface-only state variables ---
        BE%BS%n_surface = size(BE%model%surface_state_variables)

        if (BE%BS%n_surface > 0) then
            if (allocated(BE%BS%surface_state)) deallocate(BE%BS%surface_state)
            allocate(BE%BS%surface_state(BE%BS%n_surface))
            ! Allocate array for variable names
            maxlen = maxval(len_trim(BE%model%surface_state_variables%name))
            allocate(character(len=maxlen) :: BE%BS%sfcvar_names(BE%BS%n_surface))
            ! Link FABM to the surface state vector
            call BE%model%link_all_surface_state_data(BE%BS%surface_state)
            do ivar=1, BE%BS%n_surface
                BE%BS%sfcvar_names(ivar) = trim(BE%model%surface_state_variables(ivar)%name) ! Retrieve names of tracers
            end do
        end if    
        
        BE%BS%n_total = BE%BS%n_interior + BE%BS%n_surface + BE%BS%n_bottom  ! Total number of variables

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
            if (.not. allocated(BE%BS%pres)) allocate(BE%BS%pres(BE%grid%nz))
        end if
        if (BE%need_par .and. .not. allocated(BE%BS%par)) then
            allocate(BE%BS%par(BE%grid%nz))
        end if
        if (BE%need_swr .and. .not. allocated(BE%BS%swr)) then
            allocate(BE%BS%swr(BE%grid%nz))
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
        call BE%model%initialize_interior_state(1,BE%grid%nz)
        if (BE%BS%n_surface > 0) call BE%model%initialize_surface_state()
        if (BE%BS%n_bottom > 0)  call BE%model%initialize_bottom_state()

        
        BE%is_init = .true.
        write(*,'("✓ Biogeochemistry initialised successfully with ",I0," variables:")') BE%BS%n_total
        write(*,'(" ",I0," interior, ",I0," surface and ",I0," bottom variables")') BE%BS%n_interior, BE%BS%n_surface, BE%BS%n_bottom
    end subroutine init_bio_fabm

    subroutine integrate_bio_fabm()
    end subroutine integrate_bio_fabm

    subroutine end_bio_fabm()
    end subroutine end_bio_fabm


!============== INTERNAL HELPERS ===============================

    ! Provides environment data to FABM
    subroutine link_environment_data(PS, FS, BE)
        type(PhysicsState),       intent(in) :: PS
        type(ForcingSnapshot),    intent(in) :: FS
        type(BioEnv),          intent(inout) :: BE
        
        BE%BS%temp = PS%temp  !Later decide what to do with the sediment layers
        BE%BS%sal  = PS%sal
        BE%BS%rho  = PS%rho
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

end module bio_main