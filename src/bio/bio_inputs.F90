! src/bio/bio_inputs.F90
!
! Reads and stores optional biogeochemical input specifications.
! These inputs can later be used to fulfil FABM dependencies, apply
! externally prescribed source fluxes, or provide water-column relaxation targets.
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
    public :: RelaxationSpec
    public :: BioInputs

    type :: BioInputSpec
        character(:), allocatable :: key           ! YAML item name
        character(:), allocatable :: tracer        ! Name for the tracer in FABM
        character(:), allocatable :: mode          ! off/constant/file
        character(:), allocatable :: filename      ! input file
        character(:), allocatable :: name          ! variable name inside file
        character(:), allocatable :: domain        ! FABM storage domain
        character(:), allocatable :: target_domain ! Source target: water_surface/water_bottom

        real(rk) :: constant = 0.0_rk
        integer :: state_index = -1

        character(:), allocatable :: time_unit
        real(rk) :: time_scale = 1.0_rk             ! converts input flux to per second
        integer :: k_target = -1

        character(:), allocatable :: data_name
        character(:), allocatable :: time_name

        logical :: active = .false.
        logical :: has_constant = .false.
        logical :: target_found = .false.
        logical :: has_repeat_year = .false.
        integer :: repeat_year = -1
    end type BioInputSpec

    type :: RelaxationSpec
        character(:), allocatable :: key            ! YAML item name
        character(:), allocatable :: tracer         ! FABM interior tracer name
        character(:), allocatable :: mode           ! off/constant/file
        character(:), allocatable :: filename       ! input NetCDF file
        character(:), allocatable :: name           ! variable name in file
        character(:), allocatable :: depth_name     ! vertical-coordinate variable name
        character(:), allocatable :: time_name      ! time-coordinate variable name
        character(:), allocatable :: data_name      ! internal DataManager name

        real(rk) :: constant = 0.0_rk
        real(rk) :: timescale = 0.0_rk               ! user value
        character(:), allocatable :: timescale_unit
        real(rk) :: timescale_s = 0.0_rk             ! internal seconds

        integer :: state_index = -1                  ! FABM interior-state index

        logical :: active = .false.
        logical :: has_repeat_year = .false.
        integer :: repeat_year = -1
    end type RelaxationSpec

    type :: BioInputs
        logical :: is_init = .false.
        logical :: has_active_dependencies = .false.
        logical :: has_active_sources = .false.
        logical :: has_active_relaxations = .false.

        character(:), allocatable :: config_file
        type(ConfigParams) :: cfg
        type(DataManager) :: dm

        logical :: is_prepared = .false.

        real(rk), allocatable :: dep_values(:)
        real(rk), allocatable :: source_values(:)
        real(rk), allocatable :: relaxation_targets(:,:) ! (water depth, relaxation)
        integer,  allocatable :: relaxation_index(:)     ! FABM interior tracer -> relaxation entry (0 = none)

        type(BioInputSpec), allocatable :: dependencies(:)
        type(BioInputSpec), allocatable :: sources(:)
        type(RelaxationSpec), allocatable :: relaxations(:)

        character(len=PARAMLEN), allocatable :: dependency_keys(:)
        character(len=PARAMLEN), allocatable :: source_keys(:)
        character(len=PARAMLEN), allocatable :: relaxation_keys(:)

    contains
        procedure :: init         => bio_inputs_init
        procedure :: link_to_fabm => bio_inputs_link_to_fabm
        procedure :: prepare      => bio_inputs_prepare
        procedure :: tick         => bio_inputs_tick
        procedure :: update       => bio_inputs_update
        procedure :: clear        => bio_inputs_clear
    end type BioInputs

contains

    subroutine bio_inputs_init(self, input_cfg_file, FabmMod, has_input, k_wat_sfc, k_wat_btm, &
                               calendar_cfg, location, start_datetime, end_datetime, load_yearly, &
                               target_depth, ok, errmsg)
        class(BioInputs),        intent(inout) :: self
        character(*),           intent(in)    :: input_cfg_file
        class(type_fabm_model), pointer, intent(in) :: FabmMod
        logical,                intent(out)   :: has_input
        integer,                intent(in)    :: k_wat_sfc
        integer,                intent(in)    :: k_wat_btm
        type(CFCalendar),       intent(in)    :: calendar_cfg
        type(LocationInfo),     intent(in)    :: location
        type(DateTime),         intent(in)    :: start_datetime, end_datetime
        logical,                intent(in)    :: load_yearly
        real(rk), optional,     intent(in)    :: target_depth(:)
        logical, optional,      intent(out)   :: ok
        character(*), optional, intent(out)   :: errmsg

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
        call self%cfg%init()
        call self%cfg%load_yaml_content(self%config_file)

        call read_bio_input_entries(self, 'dependencies', self%dependency_keys, self%dependencies, is_source=.false.)
        call read_bio_input_entries(self, 'sources', self%source_keys, self%sources, is_source=.true.)
        call read_relaxation_entries(self, self%relaxation_keys, self%relaxations)

        self%has_active_dependencies = allocated(self%dependencies) .and. size(self%dependencies) > 0
        self%has_active_sources      = allocated(self%sources)      .and. size(self%sources) > 0
        self%has_active_relaxations  = allocated(self%relaxations)  .and. size(self%relaxations) > 0

        has_input = self%has_active_dependencies .or. self%has_active_sources .or. self%has_active_relaxations
        if (.not. has_input) return

        if (self%has_active_relaxations) then
            if (.not. present(target_depth)) then
                if (present(errmsg)) errmsg = 'BioInputs relaxation requires water-column target depths.'
                if (present(ok)) ok = .false.
                return
            end if
            if (size(target_depth) <= 0) then
                if (present(errmsg)) errmsg = 'BioInputs relaxation target depth array is empty.'
                if (present(ok)) ok = .false.
                return
            end if
        end if

        if (self%has_active_dependencies) then
            call validate_dependencies(self, FabmMod)
        end if

        if (self%has_active_sources) then
            call validate_sources(self, FabmMod, lok, msg)
            if (.not. lok) then
                call report_init_error(msg, ok, errmsg)
                return
            end if
        end if

        if (self%has_active_relaxations) then
            call validate_relaxations(self, FabmMod, lok, msg)
            if (.not. lok) then
                call report_init_error(msg, ok, errmsg)
                return
            end if
        end if

        self%has_active_dependencies = any_active_specs(self%dependencies)
        self%has_active_sources      = any_active_specs(self%sources)
        self%has_active_relaxations  = any_active_relaxations(self%relaxations)

        if (self%has_active_sources) then
            call set_target_index(self, k_wat_sfc, k_wat_btm)
        end if

        call compact_active_entries(self%dependencies)
        call compact_active_entries(self%sources)
        call compact_active_relaxations(self%relaxations)

        self%has_active_dependencies = allocated(self%dependencies) .and. size(self%dependencies) > 0
        self%has_active_sources      = allocated(self%sources)      .and. size(self%sources) > 0
        self%has_active_relaxations  = allocated(self%relaxations)  .and. size(self%relaxations) > 0

        call build_relaxation_index(self, size(FabmMod%interior_state_variables))

        has_input = self%has_active_dependencies .or. self%has_active_sources .or. self%has_active_relaxations

        if (.not. has_input) then
            write(*,'(A)') 'No active Bio input data remain after validation.'
            self%is_init = .true.
            if (present(ok)) ok = .true.
            return
        end if

        if (self%has_active_relaxations) then
            call initialise_data_manager(self, calendar_cfg, location, start_datetime, end_datetime, load_yearly, &
                                         target_depth, lok, msg)
        else
            call initialise_data_manager(self, calendar_cfg, location, start_datetime, end_datetime, load_yearly, &
                                         ok=lok, errmsg=msg)
        end if
        if (.not. lok) then
            call report_init_error(msg, ok, errmsg)
            return
        end if

        if (self%has_active_relaxations) then
            call allocate_live_storage(self, size(target_depth))
        else
            call allocate_live_storage(self)
        end if

        self%is_init = .true.
        if (present(ok)) ok = .true.
    end subroutine bio_inputs_init

    subroutine bio_inputs_prepare(self, dt_main, ok, errmsg)
        class(BioInputs), intent(inout) :: self
        integer(lk),      intent(in)    :: dt_main
        logical,          intent(out)   :: ok
        character(*),     intent(out)   :: errmsg

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
        integer(lk),      intent(in)    :: model_time
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        if (.not. self%is_prepared) return
        call self%dm%tick(model_time, ok, errmsg)
    end subroutine bio_inputs_tick

    subroutine bio_inputs_update(self, model_time, ok, errmsg)
        class(BioInputs), intent(inout) :: self
        integer(lk),      intent(in)    :: model_time
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        real(rk) :: input_value
        integer  :: i
        logical  :: lok
        character(len=512) :: lmsg

        if (present(ok)) ok = .false.
        if (present(errmsg)) errmsg = ''

        do i = 1, size(self%dependencies)
            self%dep_values(i) = self%dm%value(self%dependencies(i)%data_name, model_time, lok, lmsg)
            if (.not. lok) then
                call report_runtime_error(lmsg, ok, errmsg)
                return
            end if
        end do

        do i = 1, size(self%sources)
            input_value = self%dm%value(self%sources(i)%data_name, model_time, lok, lmsg)
            if (.not. lok) then
                call report_runtime_error(lmsg, ok, errmsg)
                return
            end if
            self%source_values(i) = input_value * self%sources(i)%time_scale
        end do

        do i = 1, size(self%relaxations)
            select case (trim(self%relaxations(i)%mode))
            case ('constant')
                input_value = self%dm%value(self%relaxations(i)%data_name, model_time, lok, lmsg)
                if (.not. lok) then
                    call report_runtime_error(lmsg, ok, errmsg)
                    return
                end if
                self%relaxation_targets(:,i) = input_value
            case ('file')
                call self%dm%profile(self%relaxations(i)%data_name, model_time, &
                                     self%relaxation_targets(:,i), lok, lmsg)
                if (.not. lok) then
                    call report_runtime_error(lmsg, ok, errmsg)
                    return
                end if
            end select
        end do

        if (present(ok)) ok = .true.
        if (present(errmsg)) errmsg = ''
    end subroutine bio_inputs_update

    subroutine bio_inputs_clear(self)
        class(BioInputs), intent(inout) :: self

        self%is_init = .false.
        self%has_active_dependencies = .false.
        self%has_active_sources = .false.
        self%has_active_relaxations = .false.
        self%is_prepared = .false.

        if (allocated(self%config_file)) deallocate(self%config_file)
        if (allocated(self%dependencies)) deallocate(self%dependencies)
        if (allocated(self%sources)) deallocate(self%sources)
        if (allocated(self%relaxations)) deallocate(self%relaxations)

        if (allocated(self%dependency_keys)) deallocate(self%dependency_keys)
        if (allocated(self%source_keys)) deallocate(self%source_keys)
        if (allocated(self%relaxation_keys)) deallocate(self%relaxation_keys)

        if (allocated(self%dep_values)) deallocate(self%dep_values)
        if (allocated(self%source_values)) deallocate(self%source_values)
        if (allocated(self%relaxation_targets)) deallocate(self%relaxation_targets)
        if (allocated(self%relaxation_index)) deallocate(self%relaxation_index)

        call self%dm%clear()
        call self%cfg%clear()
    end subroutine bio_inputs_clear

    subroutine bio_inputs_link_to_fabm(self, FabmMod)
        class(BioInputs),       intent(inout) :: self
        class(type_fabm_model), pointer, intent(in) :: FabmMod

        integer :: i
        character(:), allocatable :: name

        do i = 1, size(self%dependencies)
            name = trim(self%dependencies(i)%key)

            select case (trim(self%dependencies(i)%domain))
            case ('horizontal')
                call FabmMod%link_horizontal_data(FabmMod%get_horizontal_variable_id_by_name(name), self%dep_values(i))
            case ('scalar')
                call FabmMod%link_scalar(FabmMod%get_scalar_variable_id_by_name(name), self%dep_values(i))
            case ('interior')
                error stop 'File/constant interior dependencies are not connected yet: need depth-resolved storage.'
            end select
        end do
    end subroutine bio_inputs_link_to_fabm

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

        if (self%cfg%has_key(section)) then
            keys = self%cfg%get_child_keys(section)
        else
            allocate(keys(0))
        end if

        nactive = 0
        do i = 1, size(keys)
            base = trim(section)//'.'//trim(keys(i))
            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, &
                                          trim_value=.true., match_case=.false.)
            if (to_lower(trim(mode)) == 'file' .or. to_lower(trim(mode)) == 'constant') nactive = nactive + 1
        end do

        allocate(specs(nactive))
        nactive = 0

        do i = 1, size(keys)
            base = trim(section)//'.'//trim(keys(i))
            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, &
                                          trim_value=.true., match_case=.false.)

            select case (to_lower(trim(mode)))
            case ('off', 'false')
                cycle
            case ('file', 'constant')
                nactive = nactive + 1
                specs(nactive)%key = trim(keys(i))
                specs(nactive)%mode = to_lower(trim(mode))
                specs(nactive)%active = .true.

                if (is_source) then
                    specs(nactive)%data_name = 'src:'//trim(keys(i))
                    specs(nactive)%tracer = self%cfg%get_param_str(base//'.tracer', required=.true., trim_value=.true.)
                    specs(nactive)%target_domain = normalise_source_target( &
                        self%cfg%get_param_str(base//'.target_domain', required=.true., trim_value=.true.))
                    specs(nactive)%time_unit = normalise_time_unit( &
                        self%cfg%get_param_str(base//'.time_unit', required=.false., default='second', trim_value=.true.))
                    specs(nactive)%time_scale = 1.0_rk / seconds_per_time_unit(specs(nactive)%time_unit)
                else
                    specs(nactive)%data_name = 'dep:'//trim(keys(i))
                end if

                select case (to_lower(trim(mode)))
                case ('file')
                    specs(nactive)%filename  = self%cfg%get_param_str(base//'.filename', required=.true.)
                    specs(nactive)%name      = self%cfg%get_param_str(base//'.name', required=.true.)
                    specs(nactive)%time_name = self%cfg%get_param_str(base//'.time_name', default='time')
                    call read_optional_repeat_year(self%cfg, base, specs(nactive)%has_repeat_year, specs(nactive)%repeat_year)
                case ('constant')
                    specs(nactive)%constant = self%cfg%get_param_num(base//'.constant', finite=.true., required=.true.)
                    specs(nactive)%has_constant = .true.
                    specs(nactive)%has_repeat_year = .false.
                    specs(nactive)%repeat_year = -1
                end select
            end select
        end do
    end subroutine read_bio_input_entries

    subroutine read_relaxation_entries(self, keys, specs)
        class(BioInputs), intent(inout) :: self
        character(len=PARAMLEN), allocatable, intent(out) :: keys(:)
        type(RelaxationSpec), allocatable, intent(out) :: specs(:)

        character(len=8), dimension(4) :: mode_choices
        character(:), allocatable :: base, mode
        integer :: i, nactive

        mode_choices = ['file    ', 'constant', 'off     ', 'false   ']

        if (self%cfg%has_key('relaxation')) then
            keys = self%cfg%get_child_keys('relaxation')
        else
            allocate(keys(0))
        end if

        nactive = 0
        do i = 1, size(keys)
            base = 'relaxation.'//trim(keys(i))
            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, &
                                          trim_value=.true., match_case=.false.)
            if (to_lower(trim(mode)) == 'file' .or. to_lower(trim(mode)) == 'constant') nactive = nactive + 1
        end do

        allocate(specs(nactive))
        nactive = 0

        do i = 1, size(keys)
            base = 'relaxation.'//trim(keys(i))
            mode = self%cfg%get_param_str(base//'.mode', required=.true., choices=mode_choices, &
                                          trim_value=.true., match_case=.false.)

            select case (to_lower(trim(mode)))
            case ('off', 'false')
                cycle
            case ('file', 'constant')
                nactive = nactive + 1
                specs(nactive)%key = trim(keys(i))
                specs(nactive)%tracer = self%cfg%get_param_str(base//'.tracer', required=.true., trim_value=.true.)
                specs(nactive)%mode = to_lower(trim(mode))
                specs(nactive)%data_name = 'relax:'//trim(keys(i))
                specs(nactive)%active = .true.

                specs(nactive)%timescale = self%cfg%get_param_num(base//'.timescale', required=.true., &
                                                                  positive=.true., finite=.true.)
                specs(nactive)%timescale_unit = normalise_time_unit( &
                    self%cfg%get_param_str(base//'.timescale_unit', required=.true., trim_value=.true.))
                specs(nactive)%timescale_s = specs(nactive)%timescale * &
                                              seconds_per_time_unit(specs(nactive)%timescale_unit)

                select case (to_lower(trim(mode)))
                case ('file')
                    specs(nactive)%filename   = self%cfg%get_param_str(base//'.filename', required=.true.)
                    specs(nactive)%name       = self%cfg%get_param_str(base//'.name', required=.true.)
                    specs(nactive)%depth_name = self%cfg%get_param_str(base//'.depth_name', required=.true.)
                    specs(nactive)%time_name  = self%cfg%get_param_str(base//'.time_name', default='time')
                    call read_optional_repeat_year(self%cfg, base, specs(nactive)%has_repeat_year, specs(nactive)%repeat_year)
                case ('constant')
                    specs(nactive)%constant = self%cfg%get_param_num(base//'.constant', finite=.true., required=.true.)
                    specs(nactive)%has_repeat_year = .false.
                    specs(nactive)%repeat_year = -1
                end select
            end select
        end do
    end subroutine read_relaxation_entries

    subroutine read_optional_repeat_year(cfg, base, has_repeat_year, repeat_year)
        type(ConfigParams), intent(in) :: cfg
        character(*), intent(in) :: base
        logical, intent(out) :: has_repeat_year
        integer, intent(out) :: repeat_year

        if (.not. cfg%is_disabled(trim(base)//'.climatology_year')) then
            repeat_year = cfg%get_param_int(trim(base)//'.climatology_year')
            has_repeat_year = .true.
        else
            repeat_year = -1
            has_repeat_year = .false.
        end if
    end subroutine read_optional_repeat_year

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

    logical function any_active_relaxations(specs)
        type(RelaxationSpec), intent(in) :: specs(:)
        integer :: i
        any_active_relaxations = .false.
        do i = 1, size(specs)
            if (specs(i)%active) then
                any_active_relaxations = .true.
                return
            end if
        end do
    end function any_active_relaxations

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

            if (FabmMod%is_variable_used(FabmMod%get_horizontal_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'horizontal'
                if (FabmMod%variable_needs_values(FabmMod%get_horizontal_variable_id_by_name(name))) needs_values = .true.
            end if

            if (FabmMod%is_variable_used(FabmMod%get_interior_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'interior'
                if (FabmMod%variable_needs_values(FabmMod%get_interior_variable_id_by_name(name))) needs_values = .true.
            end if

            if (FabmMod%is_variable_used(FabmMod%get_scalar_variable_id_by_name(name))) then
                found = .true.
                self%dependencies(i)%domain = 'scalar'
                if (FabmMod%variable_needs_values(FabmMod%get_scalar_variable_id_by_name(name))) needs_values = .true.
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

            do ivar = 1, size(FabmMod%interior_state_variables)
                if (trim(name) == trim(FabmMod%interior_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'interior'
                    self%sources(i)%state_index = ivar
                end if
            end do

            do ivar = 1, size(FabmMod%surface_state_variables)
                if (trim(name) == trim(FabmMod%surface_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'surface'
                    self%sources(i)%state_index = ivar
                end if
            end do

            do ivar = 1, size(FabmMod%bottom_state_variables)
                if (trim(name) == trim(FabmMod%bottom_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    self%sources(i)%domain = 'bottom'
                    self%sources(i)%state_index = ivar
                end if
            end do

            if (nfound == 1) then
                self%sources(i)%target_found = .true.
                if (trim(self%sources(i)%domain) /= 'interior') then
                    write(*,'(A,A,A,A,A)') 'WARNING: External source "', trim(name), &
                        '" targets a FABM ', trim(self%sources(i)%domain), &
                        ' state variable. External source fluxes are currently only supported for interior variables. This source will be ignored.'
                    self%sources(i)%active = .false.
                    cycle
                end if
                self%sources(i)%active = .true.
            else if (nfound == 0) then
                errmsg = 'External input source "'//trim(self%sources(i)%key)//'" targets tracer "'//trim(name)// &
                         '", which was not found among FABM state variables.'
                return
            else
                errmsg = 'External input source "'//trim(self%sources(i)%key)//'" targets tracer "'//trim(name)// &
                         '", which matched FABM state variables in more than one domain.'
                return
            end if

            write(*,'(A,A,A,A,A,A,A,A,A)') ' - source "', trim(self%sources(i)%key), &
                '" targets tracer "', trim(name), '" which is a FABM ', trim(self%sources(i)%domain), &
                ' state variable. Source flux applied at ', trim(self%sources(i)%target_domain), '.'
        end do

        ok = .true.
    end subroutine validate_sources

    subroutine validate_relaxations(self, FabmMod, ok, errmsg)
        class(BioInputs),       intent(inout) :: self
        class(type_fabm_model), pointer, intent(in) :: FabmMod
        logical,                intent(out)   :: ok
        character(*),           intent(out)   :: errmsg

        integer :: i, j, ivar, nfound, matched_index
        character(:), allocatable :: name, domain
        logical :: disable_transport

        ok = .false.
        errmsg = ''

        if (.not. allocated(self%relaxations)) then
            ok = .true.
            return
        end if

        write(*,'(A)') 'Validating water-column relaxation inputs:'

        do i = 1, size(self%relaxations)
            name = trim(self%relaxations(i)%tracer)
            nfound = 0
            matched_index = -1
            domain = ''

            do ivar = 1, size(FabmMod%interior_state_variables)
                if (trim(name) == trim(FabmMod%interior_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    matched_index = ivar
                    domain = 'interior'
                end if
            end do

            do ivar = 1, size(FabmMod%surface_state_variables)
                if (trim(name) == trim(FabmMod%surface_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    domain = 'surface'
                end if
            end do

            do ivar = 1, size(FabmMod%bottom_state_variables)
                if (trim(name) == trim(FabmMod%bottom_state_variables(ivar)%name)) then
                    nfound = nfound + 1
                    domain = 'bottom'
                end if
            end do

            if (nfound == 0) then
                errmsg = 'Relaxation entry "'//trim(self%relaxations(i)%key)//'" targets tracer "'//trim(name)// &
                         '", which was not found among FABM state variables.'
                return
            else if (nfound > 1) then
                errmsg = 'Relaxation entry "'//trim(self%relaxations(i)%key)//'" targets tracer "'//trim(name)// &
                         '", which matched FABM state variables in more than one domain.'
                return
            else if (trim(domain) /= 'interior') then
                errmsg = 'Relaxation entry "'//trim(self%relaxations(i)%key)//'" targets tracer "'//trim(name)// &
                         '", but relaxation is supported only for FABM interior water-column state variables.'
                return
            end if

            self%relaxations(i)%state_index = matched_index

            do j = 1, i - 1
                if (self%relaxations(j)%state_index == matched_index) then
                    errmsg = 'Duplicate relaxation definitions target FABM tracer "'//trim(name)//'".'
                    return
                end if
            end do

            disable_transport = FabmMod%interior_state_variables(matched_index)%properties% &
                                get_logical('disable_transport', default=.false.)
            if (disable_transport) then
                write(*,'(A,A,A)') 'WARNING: Relaxation configured for tracer "', trim(name), &
                    '", but FABM sets disable_transport=true. Relaxation will be ignored for this tracer.'
                self%relaxations(i)%active = .false.
                cycle
            end if

            write(*,'(A,A,A,F10.3,1X,A,A)') ' - relaxation "', trim(self%relaxations(i)%key), &
                '" targets '//trim(name)//' with timescale ', self%relaxations(i)%timescale, &
                trim(self%relaxations(i)%timescale_unit), '.'
        end do

        ok = .true.
    end subroutine validate_relaxations

    function normalise_source_target(value) result(out)
        character(*), intent(in) :: value
        character(:), allocatable :: out

        select case (to_lower(trim(value)))
        case ('water_surface', 'surface', 'top')
            out = 'water_surface'
        case ('water_bottom', 'bottom')
            out = 'water_bottom'
        case default
            error stop 'Invalid source target_domain "'//trim(value)//'". Valid values are water_surface, water_bottom.'
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
            error stop 'Invalid time unit "'//trim(value)//'". Valid values are second, hour, day, year.'
        end select
    end function normalise_time_unit

    real(rk) function seconds_per_time_unit(unit) result(scale)
        character(*), intent(in) :: unit

        select case (to_lower(trim(unit)))
        case ('second')
            scale = 1.0_rk
        case ('hour')
            scale = real(sec_per_hour, rk)
        case ('day')
            scale = real(sec_per_day, rk)
        case ('year')
            scale = 365.0_rk * real(sec_per_day, rk)
        case default
            error stop 'seconds_per_time_unit: unsupported normalized time unit.'
        end select
    end function seconds_per_time_unit

    subroutine set_target_index(self, k_wat_sfc, k_wat_btm)
        class(BioInputs), intent(inout) :: self
        integer, intent(in) :: k_wat_sfc, k_wat_btm
        integer :: i

        if (.not. allocated(self%sources)) return

        do i = 1, size(self%sources)
            if (.not. self%sources(i)%active) cycle

            if (trim(self%sources(i)%domain) /= 'interior') then
                error stop 'Source "'//trim(self%sources(i)%key)//'" has unsupported FABM state domain.'
            end if

            select case (trim(self%sources(i)%target_domain))
            case ('water_surface')
                self%sources(i)%k_target = k_wat_sfc
            case ('water_bottom')
                self%sources(i)%k_target = k_wat_btm
            case default
                error stop 'Unknown source target_domain "'//trim(self%sources(i)%target_domain)//'".'
            end select
        end do
    end subroutine set_target_index

    subroutine compact_active_entries(specs)
        type(BioInputSpec), allocatable, intent(inout) :: specs(:)
        type(BioInputSpec), allocatable :: tmp(:)
        integer :: i, n, j

        if (.not. allocated(specs)) return

        n = count([(specs(i)%active, i=1,size(specs))])
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

    subroutine compact_active_relaxations(specs)
        type(RelaxationSpec), allocatable, intent(inout) :: specs(:)
        type(RelaxationSpec), allocatable :: tmp(:)
        integer :: i, n, j

        if (.not. allocated(specs)) return

        n = count([(specs(i)%active, i=1,size(specs))])
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
    end subroutine compact_active_relaxations

    subroutine build_relaxation_index(self, nint)
        class(BioInputs), intent(inout) :: self
        integer,          intent(in)    :: nint

        integer :: i, ivar

        if (allocated(self%relaxation_index)) deallocate(self%relaxation_index)
        allocate(self%relaxation_index(nint))
        self%relaxation_index = 0

        if (.not. allocated(self%relaxations)) return

        do i = 1, size(self%relaxations)
            ivar = self%relaxations(i)%state_index
            if (ivar < 1 .or. ivar > nint) then
                error stop 'build_relaxation_index: invalid FABM interior-state index.'
            end if
            self%relaxation_index(ivar) = i
        end do
    end subroutine build_relaxation_index

    subroutine initialise_data_manager(self, calendar_cfg, location, start_datetime, end_datetime, load_yearly, &
                                       target_depth, ok, errmsg)
        class(BioInputs),   intent(inout) :: self
        type(CFCalendar),   intent(in)    :: calendar_cfg
        type(LocationInfo), intent(in)    :: location
        type(DateTime),     intent(in)    :: start_datetime, end_datetime
        logical,            intent(in)    :: load_yearly
        real(rk), optional, intent(in)    :: target_depth(:)
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
        cfg%load_yearly = load_yearly

        if (self%has_active_relaxations) then
            if (.not. present(target_depth)) then
                ok = .false.
                errmsg = 'Relaxation DataSpecs require water-column target depths.'
                return
            end if
            call build_bio_data_specs(self, specs, target_depth)
        else
            call build_bio_data_specs(self, specs)
        end if

        call self%dm%init(specs, cfg, calendar_cfg, location, start_datetime, end_datetime, ok, errmsg)
    end subroutine initialise_data_manager

    subroutine build_bio_data_specs(self, specs, target_depth)
        class(BioInputs), intent(in) :: self
        type(DataSpec), allocatable, intent(out) :: specs(:)
        real(rk), optional, intent(in) :: target_depth(:)

        integer :: n, i, j

        n = size(self%dependencies) + size(self%sources) + size(self%relaxations)
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

        do i = 1, size(self%relaxations)
            j = j + 1
            if (.not. present(target_depth)) error stop 'build_bio_data_specs: missing target depths for relaxation.'
            call spec_from_relaxation_entry(self%relaxations(i), target_depth, specs(j))
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
            spec%path = trim(entry%filename)
            spec%time_var = trim(entry%time_name)
            call apply_repeat_metadata(entry%has_repeat_year, entry%repeat_year, spec)
        case ('constant')
            spec%input_type = DATA_INPUT_CONSTANT
            spec%const_value = entry%constant
            spec%source_var = ''
            spec%path = ''
            spec%time_var = ''
        end select
    end subroutine spec_from_bio_entry

    subroutine spec_from_relaxation_entry(entry, target_depth, spec)
        type(RelaxationSpec), intent(in) :: entry
        real(rk), intent(in) :: target_depth(:)
        type(DataSpec), intent(out) :: spec

        spec%name = trim(entry%data_name)

        select case (trim(entry%mode))
        case ('file')
            spec%input_type = DATA_INPUT_FILE
            spec%source_var = trim(entry%name)
            spec%path = trim(entry%filename)
            spec%time_var = trim(entry%time_name)
            spec%is_profile = .true.
            spec%depth_var = trim(entry%depth_name)
            spec%target_depth = target_depth
            call apply_repeat_metadata(entry%has_repeat_year, entry%repeat_year, spec)
        case ('constant')
            spec%input_type = DATA_INPUT_CONSTANT
            spec%const_value = entry%constant
            spec%source_var = ''
            spec%path = ''
            spec%time_var = ''
            spec%is_profile = .false.
        end select
    end subroutine spec_from_relaxation_entry

    subroutine apply_repeat_metadata(has_repeat_year, repeat_year, spec)
        logical, intent(in) :: has_repeat_year
        integer, intent(in) :: repeat_year
        type(DataSpec), intent(inout) :: spec

        if (has_repeat_year) then
            spec%repeat_enabled = .true.
            spec%repeat_year = repeat_year
            spec%time_mode = DATA_TIME_REPEAT_YEAR
        else
            spec%repeat_enabled = .false.
            spec%repeat_year = -huge(1)
            spec%time_mode = DATA_TIME_ABSOLUTE
        end if
    end subroutine apply_repeat_metadata

    subroutine allocate_live_storage(self, nwat)
        class(BioInputs), intent(inout) :: self
        integer, optional, intent(in) :: nwat

        if (allocated(self%dep_values)) deallocate(self%dep_values)
        if (allocated(self%source_values)) deallocate(self%source_values)
        if (allocated(self%relaxation_targets)) deallocate(self%relaxation_targets)

        allocate(self%dep_values(size(self%dependencies)))
        allocate(self%source_values(size(self%sources)))
        self%dep_values = 0.0_rk
        self%source_values = 0.0_rk

        if (size(self%relaxations) > 0) then
            if (.not. present(nwat)) error stop 'allocate_live_storage: missing water-column size for relaxation.'
            allocate(self%relaxation_targets(nwat, size(self%relaxations)))
        else
            allocate(self%relaxation_targets(0, 0))
        end if
        self%relaxation_targets = 0.0_rk
    end subroutine allocate_live_storage

    subroutine report_init_error(msg, ok, errmsg)
        character(*), intent(in) :: msg
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        if (present(errmsg)) errmsg = trim(msg)
        if (present(ok)) then
            ok = .false.
        else
            write(*,*) trim(msg)
            stop 1
        end if
    end subroutine report_init_error

    subroutine report_runtime_error(msg, ok, errmsg)
        character(*), intent(in) :: msg
        logical, optional, intent(out) :: ok
        character(*), optional, intent(out) :: errmsg

        if (present(errmsg)) then
            errmsg = trim(msg)
        else
            write(*,*) 'ERROR: ', trim(msg)
            stop 1
        end if
        if (present(ok)) ok = .false.
    end subroutine report_runtime_error

end module bio_inputs