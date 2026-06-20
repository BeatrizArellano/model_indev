! src/bio/bio_inputs.F90
!
! Reads and stores optional biogeochemical input specifications.
! These inputs can later be used to fulfil FABM dependencies or apply
! externally prescribed source fluxes.
!
module bio_inputs
    use data_manager,     only: DataManager
    use data_types,       only: DataLoaderCfg, DataSpec, &
                                DATA_INPUT_FILE, DATA_INPUT_CONSTANT, &
                                DATA_TIME_ABSOLUTE, DATA_TIME_REPEAT_YEAR
    use fabm,             only: type_fabm_model
    use geo_utils,        only: LocationInfo
    use precision_types,  only: rk, lk
    use read_config_yaml, only: ConfigParams, PARAMLEN
    use str_utils,        only: to_lower  
    use time_types,       only: DateTime, CFCalendar, cal_unknown
    use time_utils,       only: sec_per_day, sec_per_hour

    implicit none
    private

    public :: BioInputSpec
    public :: BioInputs

    type :: BioInputSpec 
        character(:), allocatable :: key           ! YAML item name 
        character(:), allocatable :: tracer        ! Name for the tracer in FABM
        character(:), allocatable :: mode          ! off/constant/file
        character(:), allocatable :: filename      ! input file
        character(:), allocatable :: name          ! variable name inside file
        character(:), allocatable :: domain        ! FABM storage domain: interior/surface/bottom/horizontal/scalar
        character(:), allocatable :: target_domain ! Source application target: water_surface/water_bottom

        real(rk) :: constant = 0.0_rk
        ! FABM index for state variables (Only applicable to sources)
        integer :: state_index = -1   
        
        character(:), allocatable :: time_unit
        real(rk) :: time_scale = 1.0_rk   ! converts input flux to per second
        integer :: k_target = -1   ! Index in the vertical grid where the source flux will be applied

        character(:), allocatable :: data_name
        character(:), allocatable :: time_name

        logical :: active = .false.   
        logical :: has_constant = .false.
        logical :: target_found = .false.
        logical :: has_repeat_year = .false.
        integer :: repeat_year = -1
    contains
       !procedure :: clear => clear_bio_input_spec
    end type BioInputSpec
 
 
    type :: BioInputs
        logical :: is_init = .false.
        logical :: has_active_dependencies = .false.
        logical :: has_active_sources = .false.

        character(:), allocatable :: config_file
        type(ConfigParams) :: cfg

        type(DataManager) :: dm

        logical :: is_prepared = .false.

        real(rk), allocatable :: dep_values(:)
        real(rk), allocatable :: source_values(:)
    
        type(BioInputSpec), allocatable :: dependencies(:)
        type(BioInputSpec), allocatable :: sources(:)

        character(len=PARAMLEN), allocatable  :: dependency_keys(:)
        character(len=PARAMLEN), allocatable :: source_keys(:)
 
    contains
        procedure :: init         => bio_inputs_init
        procedure :: link_to_fabm => bio_inputs_link_to_fabm
        procedure :: prepare      => bio_inputs_prepare   ! set up the DataManager for the simulation timestep
        procedure :: tick         => bio_inputs_tick      ! let the DataManager load/advance data for the current model time
        procedure :: update       => bio_inputs_update    ! copy sampled values into BioInputs arrays
        procedure :: clear        => bio_inputs_clear
    end type BioInputs

contains

    subroutine bio_inputs_init(self, input_cfg_file, FabmMod, has_input, k_wat_sfc, k_wat_btm, &                               
                               calendar_cfg, location, start_datetime, end_datetime, load_yearly, ok, errmsg)
        class(BioInputs),        intent(inout) :: self
        character(*),            intent(in)    :: input_cfg_file
        class(type_fabm_model),  pointer, intent(in) :: FabmMod
        logical,                 intent(out)   :: has_input
        integer,                 intent(in)    :: k_wat_sfc
        integer,                 intent(in)    :: k_wat_btm
        type(CFCalendar),        intent(in)    :: calendar_cfg
        type(LocationInfo),      intent(in)    :: location
        type(DateTime),          intent(in)    :: start_datetime, end_datetime
        logical,                 intent(in)    :: load_yearly
        logical,       optional, intent(out)   :: ok
        character(*),  optional, intent(out)   :: errmsg

        logical :: lok
        character(len=512) :: msg

        lok = .false.    
        
        has_input = .false.
        if (present(ok)) ok = .false.

        if (present(errmsg)) errmsg = ''

        call self%clear()

        if (.not. associated(FabmMod)) then
            if (present(errmsg)) errmsg = 'FABM model pointer is not associated.'
            if (present(ok)) ok = .false.
            return
        end if

        if (len_trim(input_cfg_file) == 0) then
            if (present(errmsg)) errmsg = 'Bio input configuration file path is empty.'
            if (present(ok)) ok = .false.
            return
        end if

        self%config_file = trim(input_cfg_file)
        ! Read configuration file
        call self%cfg%init()
        call self%cfg%load_yaml_content(self%config_file)

        ! -------------------------------
        ! Read dependency entries
        ! -------------------------------
        call read_bio_input_entries(self, 'dependencies', self%dependency_keys, self%dependencies, is_source=.false.)
        
        ! -------------------------------
        ! Read source entries
        ! -------------------------------
        call read_bio_input_entries(self, 'sources', self%source_keys, self%sources, is_source=.true.)

        self%has_active_dependencies = allocated(self%dependencies) .and. size(self%dependencies) > 0
        self%has_active_sources      = allocated(self%sources)      .and. size(self%sources) > 0

        has_input = self%has_active_dependencies .or. self%has_active_sources

        if (.not. has_input) return

        !---------------------------------------
        !  Validate dependencies with FABM
        !---------------------------------------
        if (self%has_active_dependencies) then
            call validate_dependencies(self, FabmMod)
        end if

        !---------------------------------------
        !  Validate sources against FABM state variables
        !---------------------------------------
        if (self%has_active_sources) then
            call validate_sources(self, FabmMod, lok, msg)
            if (.not. lok) then
                if (present(errmsg)) errmsg = trim(msg)
                if (present(ok)) then
                    ok = .false.
                    return
                else
                    write(*,*) trim(msg)
                    stop 1
                end if
            end if
        end if
        ! Check again if there are active sources or dependencies
        self%has_active_dependencies = any_active_specs(self%dependencies)
        self%has_active_sources      = any_active_specs(self%sources)   

        if (self%has_active_sources) then
            call set_target_index(self, k_wat_sfc, k_wat_btm)
        end if
        
        ! Keep only entries that survived validation
        call compact_active_entries(self%dependencies)
        call compact_active_entries(self%sources)

        self%has_active_dependencies = allocated(self%dependencies) .and. size(self%dependencies) > 0
        self%has_active_sources      = allocated(self%sources)      .and. size(self%sources) > 0

        has_input = self%has_active_dependencies .or. self%has_active_sources

        if (.not. has_input) then
            write(*,'(A)') 'No active Bio input data remain after validation.'
            self%is_init = .true.
            if (present(ok)) ok = .true.
            return
        end if

        call initialise_data_manager(self, calendar_cfg, location, start_datetime, end_datetime, load_yearly, lok, msg)
        if (.not. lok) then
            if (present(errmsg)) errmsg = trim(msg)
            if (present(ok)) ok = .false.
            return
        end if

        call allocate_live_storage(self)

        self%is_init = .true.
        lok = .true.
        if (present(ok)) ok = lok        
    end subroutine bio_inputs_init

    subroutine bio_inputs_prepare(self, dt_main, ok, errmsg)
        class(BioInputs), intent(inout) :: self
        integer(lk), intent(in) :: dt_main
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        if (.not. self%is_init) then
            ok = .false.
            errmsg = 'BioInputs not initialized.'
            return
        end if

        call self%dm%prepare(dt_main, ok, errmsg)
        if (ok) self%is_prepared = .true.
    end subroutine bio_inputs_prepare

    subroutine bio_inputs_tick(self, model_time, ok, errmsg)
        class(BioInputs), intent(inout) :: self
        integer(lk), intent(in) :: model_time
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        if (.not. self%is_prepared) return

        call self%dm%tick(model_time, ok, errmsg)
    end subroutine bio_inputs_tick

    subroutine bio_inputs_update(self, model_time, ok, errmsg)
        class(BioInputs), intent(inout) :: self
        integer(lk), intent(in) :: model_time
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        real(rk) :: input_value
        integer  :: i
        logical  :: lok
        character(len=512) :: lmsg

        if (present(ok)) ok = .false.
        if (present(errmsg)) errmsg = ''

        do i = 1, size(self%dependencies)

            self%dep_values(i) = self%dm%value( &
                                    self%dependencies(i)%data_name, &
                                    model_time, lok, lmsg)
            if (.not. lok) then
                if (present(errmsg)) then
                    errmsg = trim(lmsg)
                else
                    write(*,*) 'ERROR: ', trim(lmsg)
                    stop 1
                end if
                return
            end if
        end do

        do i = 1, size(self%sources)
            input_value = self%dm%value(self%sources(i)%data_name, model_time, lok, lmsg)

            if (.not. lok) then
                if (present(errmsg)) then
                    errmsg = trim(lmsg)
                else
                    write(*,*) 'ERROR: ', trim(lmsg)
                    stop 1
                end if
                return
            end if

            ! Multiply by the time_scale to convert to flux per second
            self%source_values(i) = input_value * self%sources(i)%time_scale
        end do

        if (present(ok)) ok = .true.
        if (present(errmsg)) errmsg = ''

    end subroutine bio_inputs_update


    subroutine bio_inputs_clear(self)
        class(BioInputs), intent(inout) :: self

        self%is_init = .false.
        self%has_active_dependencies = .false.
        self%has_active_sources = .false.
        self%is_prepared = .false.

        if (allocated(self%config_file))  deallocate(self%config_file)
        if (allocated(self%dependencies)) deallocate(self%dependencies)
        if (allocated(self%sources))      deallocate(self%sources)

        if (allocated(self%dependency_keys)) deallocate(self%dependency_keys)
        if (allocated(self%source_keys))     deallocate(self%source_keys)

        if (allocated(self%dep_values)) deallocate(self%dep_values)
        if (allocated(self%source_values)) deallocate(self%source_values)

        call self%dm%clear()
        call self%cfg%clear()

    end subroutine bio_inputs_clear

    ! Link dependencies to FABM
    subroutine bio_inputs_link_to_fabm(self, FabmMod)
        class(BioInputs), intent(inout) :: self
        class(type_fabm_model), pointer, intent(in) :: FabmMod

        integer :: i
        character(:), allocatable :: name

        do i = 1, size(self%dependencies)
            name = trim(self%dependencies(i)%key)

            select case (trim(self%dependencies(i)%domain))
            case ('horizontal')
                call FabmMod%link_horizontal_data(FabmMod%get_horizontal_variable_id_by_name(name), &
                                                  self%dep_values(i))

            case ('scalar')
                call FabmMod%link_scalar(FabmMod%get_scalar_variable_id_by_name(name), &
                                         self%dep_values(i))

            case ('interior')
                error stop 'File/constant interior dependencies are not connected yet: need depth-resolved storage.'
            end select
        end do
    end subroutine bio_inputs_link_to_fabm


    !===========================================================================
    !                               Local routines
    !===========================================================================

    subroutine read_bio_input_entries(self, section, keys, specs, is_source)
        class(BioInputs), intent(inout) :: self
        character(*), intent(in) :: section
        character(len=PARAMLEN), allocatable, intent(out) :: keys(:)
        type(BioInputSpec), allocatable, intent(out) :: specs(:)
        logical, intent(in) :: is_source

        character(len=8), dimension(4) :: mode_choices
        character(:), allocatable :: base, mode
        integer :: i, nactive
        

        mode_choices = ['file    ', 'constant', 'off     ', 'false   ']

        ! -------------------------------
        ! List input entries
        ! -------------------------------
        if (self%cfg%has_key(section)) then
            keys = self%cfg%get_child_keys(section)
        else
            allocate(keys(0))
        end if

        ! -------------------------------------
        ! List active variables (mode not off)
        ! -------------------------------------
        nactive = 0
        do i = 1, size(keys)
            base = trim(section)//'.'//trim(keys(i))

            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, trim_value=.true., match_case=.false.)                   

            select case (to_lower(trim(mode)))
                case ('file', 'constant')
                    nactive = nactive + 1
                case ('off', 'false')
                    cycle
            end select
        end do

        allocate(specs(nactive))

        ! -------------------------------
        ! Fill in the InputSpec type
        ! -------------------------------
        nactive = 0
        do i = 1, size(keys)
            base = trim(section)//'.'//trim(keys(i))

            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, trim_value=.true., match_case=.false.)                

            select case (to_lower(trim(mode)))

            case ('off', 'false')
                cycle

            case ('file', 'constant')
                nactive = nactive + 1

                specs(nactive)%key    = trim(keys(i))
                specs(nactive)%mode   = to_lower(trim(mode))
                specs(nactive)%active = .true.

                if (is_source) then
                    specs(nactive)%data_name = 'src:'//trim(keys(i))
                else
                    specs(nactive)%data_name = 'dep:'//trim(keys(i))
                end if

                if (is_source) then
                    specs(nactive)%tracer = self%cfg%get_param_str(base//'.tracer', required=.true., trim_value=.true.)
                    specs(nactive)%target_domain = normalise_source_target(self%cfg%get_param_str(base//'.target_domain', required=.true., trim_value=.true.))
                    specs(nactive)%time_unit = normalise_time_unit(self%cfg%get_param_str(base//'.time_unit', required=.false., default='second', trim_value=.true.))
                    select case (to_lower(trim(specs(nactive)%time_unit)))
                        case ('second')
                            specs(nactive)%time_scale = 1.0_rk
                        case ('hour')
                            specs(nactive)%time_scale = 1.0_rk / real(sec_per_hour, rk)
                        case ('day')
                            specs(nactive)%time_scale = 1.0_rk / real(sec_per_day, rk)
                        case ('year')
                            specs(nactive)%time_scale = 1.0_rk / (365.0_rk * real(sec_per_day, rk))
                    end select
                end if

                select case (to_lower(trim(mode)))
                    case ('file')
                        specs(nactive)%filename  = self%cfg%get_param_str(base//'.filename', required=.true.)
                        specs(nactive)%name      = self%cfg%get_param_str(base//'.name',     required=.true.)
                        specs(nactive)%time_name = self%cfg%get_param_str(base//'.time_name', default='time')
                        if (.not. self%cfg%is_disabled(base//'.climatology_year')) then
                            specs(nactive)%repeat_year = self%cfg%get_param_int(base//'.climatology_year')
                            specs(nactive)%has_repeat_year = .true.
                        else
                            specs(nactive)%repeat_year = -1
                            specs(nactive)%has_repeat_year = .false.
                        end if
                    case ('constant')
                        specs(nactive)%constant = self%cfg%get_param_num(base//'.constant', finite=.true., required=.true.)
                        specs(nactive)%has_constant = .true.
                        specs(nactive)%has_repeat_year = .false.
                        specs(nactive)%repeat_year = -1
                end select
            end select
        end do
    end subroutine read_bio_input_entries

    logical function any_active_specs(specs)
        type(BioInputSpec), intent(in) :: specs(:)
        integer :: i

        any_active_specs = .false.

        do i = 1, size(specs)
            if (specs(i)%active) then
                any_active_specs = .true.
                return
            end if
        end do
    end function any_active_specs

    subroutine validate_dependencies(self, FabmMod)
        class(BioInputs),       intent(inout) :: self
        class(type_fabm_model), pointer, intent(in) :: FabmMod

        integer :: i
        logical :: found, needs_values
        character(:), allocatable :: name

        if (.not. allocated(self%dependencies)) return
        if (size(self%dependencies) == 0) return

        write(*,'(A)') 'Validating dependencies set in the input configuration file against FABM:'

        do i = 1, size(self%dependencies)
            name = trim(self%dependencies(i)%key)

            found = .false.
            needs_values = .false.

            ! -------------------------------
            ! Horizontal dependency
            ! -------------------------------
            if (FabmMod%is_variable_used(FabmMod%get_horizontal_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'horizontal'

                if (FabmMod%variable_needs_values(FabmMod%get_horizontal_variable_id_by_name(name))) then
                    needs_values = .true.
                end if
            end if

            ! -------------------------------
            ! Interior dependency
            ! -------------------------------
            if (FabmMod%is_variable_used(FabmMod%get_interior_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'interior'

                if (FabmMod%variable_needs_values(FabmMod%get_interior_variable_id_by_name(name))) then
                    needs_values = .true.
                end if
            end if

            ! -------------------------------
            ! Scalar dependency
            ! -------------------------------
            if (FabmMod%is_variable_used(FabmMod%get_scalar_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'scalar'

                if (FabmMod%variable_needs_values(FabmMod%get_scalar_variable_id_by_name(name))) then
                    needs_values = .true.
                end if
            end if

            if (.not. found) then
                write(*,'(A,A,A)') 'WARNING: Dependency "', trim(name), &
                    '" was configured but was not found as an active FABM dependency.'
                self%dependencies(i)%active = .false.
            else if (.not. needs_values) then
                write(*,'(A,A,A)') 'WARNING: Dependency "', trim(name), &
                    '" was found in FABM but does not require externally supplied values. No values will be supplied for this variable.'
                self%dependencies(i)%active = .false.

            else if (trim(self%dependencies(i)%domain) == 'interior') then
                write(*,'(A,A,A)') 'WARNING: Interior dependency "', trim(name), &
                    '" is configured but BioInputs does not yet support depth-resolved interior input. This dependency will not be supplied here.'
                self%dependencies(i)%active = .false.
            else
                write(*,'(A,A,A,A)') ' - ', trim(name), ' is a FABM ', trim(self%dependencies(i)%domain)//' dependency.'
            end if

        end do
    end subroutine validate_dependencies

    subroutine validate_sources(self, FabmMod, ok, errmsg)
        class(BioInputs),       intent(inout) :: self
        class(type_fabm_model), pointer, intent(in) :: FabmMod
        logical,                intent(out)   :: ok
        character(*),           intent(out)   :: errmsg


        integer :: i, ivar, nfound
        character(:), allocatable :: name

        errmsg = ''
        ok = .false.

        if (.not. allocated(self%sources)) then
            ok = .true.
            return
        end if    

        write(*,'(A)') 'Validating external source flux inputs:'

        do i = 1, size(self%sources)
            name = trim(self%sources(i)%tracer)
            nfound = 0

            ! -------------------------------
            ! Interior state variables
            ! -------------------------------
            do ivar = 1, size(FabmMod%interior_state_variables)
                if (trim(name) == trim(FabmMod%interior_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'interior'
                    self%sources(i)%state_index  = ivar
                end if
            end do

            ! -------------------------------
            ! Surface state variables
            ! -------------------------------
            do ivar = 1, size(FabmMod%surface_state_variables)
                if (trim(name) == trim(FabmMod%surface_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'surface'
                    self%sources(i)%state_index  = ivar
                end if
            end do

            ! -------------------------------
            ! Bottom state variables
            ! -------------------------------
            do ivar = 1, size(FabmMod%bottom_state_variables)
                if (trim(name) == trim(FabmMod%bottom_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'bottom'
                    self%sources(i)%state_index  = ivar
                end if
            end do

            if (nfound == 1) then
                self%sources(i)%target_found = .true.
                ! --------------------------------------------------------
                ! External source fluxes are currently only supported
                ! for interior FABM state variables.
                !
                ! Surface and bottom FABM variables are area-attached
                ! quantities and might not treated as volumetric tracers.
                ! --------------------------------------------------------
                if (trim(self%sources(i)%domain) /= 'interior') then
                    write(*,'(A,A,A,A,A)') &
                        'WARNING: External source "', trim(name), &
                        '" targets a FABM ', trim(self%sources(i)%domain), &
                        ' state variable. External source fluxes are currently only supported for interior variables. This source will be ignored.'
                    self%sources(i)%active = .false.
                    cycle
                end if
                self%sources(i)%active = .true.

            else if (nfound == 0) then
                errmsg = 'External input source "'//trim(self%sources(i)%key)// &
                         '" targets tracer "'//trim(name)// &
                         '", which was not found among FABM state variables.'
                return

            else if (nfound > 1) then
                errmsg = 'External input source "'//trim(self%sources(i)%key)// &
                        '" targets tracer "'//trim(name)// &
                        '", which matched FABM state variables in more than one domain.'
                return
            end if

            write(*,'(A,A,A,A,A,A,A,A,A)') ' - source "', trim(self%sources(i)%key), &
                    '" targets tracer "', trim(name), '" which is a FABM ', &
                    trim(self%sources(i)%domain), ' state variable. Source flux applied at ', &
                    trim(self%sources(i)%target_domain), '.'                   
            
        end do

        ok = .true.
    end subroutine validate_sources

    function normalise_source_target(value) result(out)
        character(*), intent(in) :: value
        character(:), allocatable :: out

        select case (to_lower(trim(value)))
            case ('water_surface', 'surface', 'top')
                out = 'water_surface'

            case ('water_bottom', 'bottom')
                out = 'water_bottom'

            case default
                error stop 'Invalid source target_domain "'//trim(value)// &
                    '". Valid values are water_surface, water_bottom.'
        end select
    end function normalise_source_target

    function normalise_time_unit(value) result(out)
        character(*), intent(in) :: value
        character(:), allocatable :: out

        select case (to_lower(trim(value)))
        case ('second', 'seconds', 'sec', 'secs', 's')
            out = 'second'

        case ('hour', 'hours', 'hr', 'hrs', 'h')
            out = 'hour'

        case ('day', 'days', 'd')
            out = 'day'

        case ('year', 'years', 'yr', 'yrs', 'y')
            out = 'year'

        case default
            error stop 'Invalid source time_unit "'//trim(value)// &
                '". Valid values are second, hour, day, year.'
        end select
    end function normalise_time_unit

    subroutine set_target_index(self, k_wat_sfc, k_wat_btm)
        class(BioInputs), intent(inout) :: self
        integer, intent(in) :: k_wat_sfc, k_wat_btm
        integer :: i

        if (.not. allocated(self%sources)) return

        do i = 1, size(self%sources)
            if (.not. self%sources(i)%active) cycle

            select case (trim(self%sources(i)%domain))

                case ('interior')
                    ! Interior variables can receive fluxes at water_surface or water_bottom.
                case default
                    error stop 'Source "'//trim(self%sources(i)%key)// &
                        '" has unsupported FABM state domain "'//trim(self%sources(i)%domain)//'".'
            end select

            select case (trim(self%sources(i)%target_domain))

                case ('water_surface')
                    self%sources(i)%k_target = k_wat_sfc

                case ('water_bottom')
                    self%sources(i)%k_target = k_wat_btm

                case default
                    error stop 'Unknown source target_domain "'// &
                        trim(self%sources(i)%target_domain)//'".'
            end select
        end do
    end subroutine set_target_index

    subroutine compact_active_entries(specs)
        type(BioInputSpec), allocatable, intent(inout) :: specs(:)

        type(BioInputSpec), allocatable :: tmp(:)
        integer :: i, n, j

        if (.not. allocated(specs)) return

        n = 0
        do i = 1, size(specs)
            if (specs(i)%active) n = n + 1
        end do

        if (n == 0) then
            deallocate(specs)
            allocate(specs(0))
            return
        end if

        allocate(tmp(n))

        j = 0
        do i = 1, size(specs)
            if (.not. specs(i)%active) cycle
            j = j + 1
            tmp(j) = specs(i)
        end do

        call move_alloc(tmp, specs)
    end subroutine compact_active_entries

    subroutine initialise_data_manager(self, calendar_cfg, location, start_datetime, end_datetime, load_yearly, ok, errmsg)
        class(BioInputs),   intent(inout) :: self
        type(CFCalendar),   intent(in)    :: calendar_cfg
        type(LocationInfo), intent(in)    :: location
        type(DateTime),     intent(in)    :: start_datetime, end_datetime
        logical,            intent(in)    :: load_yearly
        logical,            intent(out)   :: ok
        character(*),       intent(out)   :: errmsg

        type(DataLoaderCfg) :: cfg
        type(DataSpec), allocatable :: specs(:)

        if (calendar_cfg%kind == cal_unknown) then
            ok = .false.
            errmsg = 'BioInputs requires a known simulation calendar; it cannot derive the calendar from biogeochemical input files.'
            return
        end if

        cfg%cfg_calendar = calendar_cfg%kind
        cfg%load_yearly  = load_yearly

        call build_bio_data_specs(self, specs)

        call self%dm%init(specs, cfg, calendar_cfg, location, start_datetime, end_datetime, ok, errmsg)
    end subroutine initialise_data_manager

    subroutine build_bio_data_specs(self, specs)
        class(BioInputs), intent(in) :: self
        type(DataSpec), allocatable, intent(out) :: specs(:)

        integer :: n, i, j

        n = size(self%dependencies) + size(self%sources)
        allocate(specs(n))

        j = 0

        do i = 1, size(self%dependencies)
            j = j + 1
            call spec_from_bio_entry(self%dependencies(i), specs(j))
        end do

        do i = 1, size(self%sources)
            j = j + 1
            call spec_from_bio_entry(self%sources(i), specs(j))
        end do
    end subroutine build_bio_data_specs

    subroutine spec_from_bio_entry(entry, spec)
        type(BioInputSpec), intent(in) :: entry
        type(DataSpec), intent(out) :: spec

        spec%name = trim(entry%data_name)

        select case (trim(entry%mode))
        case ('file')
            spec%input_type = DATA_INPUT_FILE
            spec%source_var = trim(entry%name)
            spec%path       = trim(entry%filename)
            spec%time_var   = trim(entry%time_name)

            if (entry%has_repeat_year) then
                spec%repeat_enabled = .true.
                spec%repeat_year    = entry%repeat_year
                spec%time_mode      = DATA_TIME_REPEAT_YEAR
            else
                spec%repeat_enabled = .false.
                spec%repeat_year    = -huge(1)
                spec%time_mode      = DATA_TIME_ABSOLUTE
            end if

        case ('constant')
            spec%input_type  = DATA_INPUT_CONSTANT
            spec%const_value = entry%constant
            spec%source_var  = ''
            spec%path        = ''
            spec%time_var    = ''
        end select
    end subroutine spec_from_bio_entry

    subroutine allocate_live_storage(self)
        class(BioInputs), intent(inout) :: self

        if (allocated(self%dep_values)) deallocate(self%dep_values)
        if (allocated(self%source_values)) deallocate(self%source_values)

        allocate(self%dep_values(size(self%dependencies)))
        allocate(self%source_values(size(self%sources)))

        self%dep_values = 0.0_rk
        self%source_values = 0.0_rk
    end subroutine allocate_live_storage

end module bio_inputs

