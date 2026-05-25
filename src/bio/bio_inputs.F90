! src/bio/bio_inputs.F90
!
! Reads and stores optional biogeochemical input specifications.
! These inputs can later be used to fulfil FABM dependencies or apply
! externally prescribed source fluxes.
!
module bio_inputs
    use fabm,             only: type_fabm_model
    use precision_types,  only: rk
    use read_config_yaml, only: ConfigParams, PARAMLEN
    use str_utils,        only: to_lower  
    use time_utils,       only: sec_per_day, sec_per_hour

    implicit none
    private

    public :: BioInputSpec
    public :: BioInputs

    type :: BioInputSpec 
        character(:), allocatable :: key           ! YAML item name / FABM-facing name
        character(:), allocatable :: mode          ! off/constant/file
        character(:), allocatable :: filename      ! input file
        character(:), allocatable :: name          ! variable name inside file
        character(:), allocatable :: domain        ! FABM storage domain: interior/surface/bottom/horizontal/scalar
        character(:), allocatable :: target_domain ! Source application target: water_surface/water_bottom/sediment_surface

        real(rk) :: constant = 0.0_rk
        ! FABM index for state variables (Only applicable to sources)
        integer :: state_index = -1   
        
        character(:), allocatable :: time_unit
        real(rk) :: time_scale = 1.0_rk   ! converts input flux to per second
        integer :: k_target = -1   ! Index in the vertical grid where the source flux will be applied

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
 
       type(BioInputSpec), allocatable :: dependencies(:)
       type(BioInputSpec), allocatable :: sources(:)

       character(len=PARAMLEN), allocatable  :: dependency_keys(:)
       character(len=PARAMLEN), allocatable :: source_keys(:)
 
    contains
       procedure :: init  => bio_inputs_init
       procedure :: clear => bio_inputs_clear
    end type BioInputs

contains

    subroutine bio_inputs_init(self, input_cfg_file, FabmMod, has_input, &
                               sediments_enabled, k_wat_sfc, k_wat_btm, k_sed_sfc, &
                               ok, errmsg)
        class(BioInputs),        intent(inout) :: self
        character(*),            intent(in)    :: input_cfg_file
        class(type_fabm_model),  pointer, intent(in) :: FabmMod
        logical,                 intent(out)   :: has_input
        logical,                 intent(in)    :: sediments_enabled
        integer,                 intent(in)    :: k_wat_sfc
        integer,                 intent(in)    :: k_wat_btm
        integer,                 intent(in)    :: k_sed_sfc
        logical,       optional, intent(out)   :: ok
        character(*),  optional, intent(out)   :: errmsg

        logical :: lok
        character(len=512) :: msg

        lok = .false.        

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
            call set_target_index(self, sediments_enabled, k_wat_sfc, k_wat_btm, k_sed_sfc)
        end if
        
        ! Keep only entries that survived validation
        call compact_active_entries(self%dependencies)
        call compact_active_entries(self%sources)

        self%has_active_dependencies = allocated(self%dependencies) .and. size(self%dependencies) > 0
        self%has_active_sources      = allocated(self%sources)      .and. size(self%sources) > 0

        has_input = self%has_active_dependencies .or. self%has_active_sources

        self%is_init = .true.
        lok = .true.
        if (present(ok)) ok = lok

        if (.not. has_input) then
            write(*,'(A)') 'No active Bio input data remain after validation.'
        end if
    end subroutine bio_inputs_init


    subroutine bio_inputs_clear(self)
        class(BioInputs), intent(inout) :: self

        self%is_init = .false.
        self%has_active_dependencies = .false.
        self%has_active_sources = .false.

        if (allocated(self%config_file))  deallocate(self%config_file)
        if (allocated(self%dependencies)) deallocate(self%dependencies)
        if (allocated(self%sources))      deallocate(self%sources)

        if (allocated(self%dependency_keys)) deallocate(self%dependency_keys)
        if (allocated(self%source_keys))     deallocate(self%source_keys)

        call self%cfg%clear()

    end subroutine bio_inputs_clear


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
                        specs(nactive)%filename = self%cfg%get_param_str(base//'.filename', required=.true.)
                        specs(nactive)%name     = self%cfg%get_param_str(base//'.name',     required=.true.)
                    case ('constant')
                        specs(nactive)%constant = self%cfg%get_param_num(base//'.constant', finite=.true., required=.true.)
                        specs(nactive)%has_constant = .true.
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
            name = trim(self%sources(i)%key)
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
                self%sources(i)%active = .true.

            else if (nfound == 0) then
                errmsg = 'External input source "'//trim(name)//'" was not found among FABM state variables.'
                return

            else if (nfound > 1) then
                errmsg = 'External input source "'//trim(name)//'" matched FABM state variables in more than one domain.'
                return
            end if

            write(*,'(A,A,A,A,A,A,A)') ' - ', trim(name), ' is a FABM ', trim(self%sources(i)%domain), &
                  ' state variable. Source flux to be applied at the ', &
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

            case ('sediment_surface', 'sed_surface', 'sediment_top')
                out = 'sediment_surface'

            case default
                error stop 'Invalid source target_domain "'//trim(value)// &
                    '". Valid values are water_surface, water_bottom, sediment_surface.'
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

    subroutine set_target_index(self, sediments_enabled, k_wat_sfc, k_wat_btm, k_sed_sfc)
        class(BioInputs), intent(inout) :: self
        logical, intent(in) :: sediments_enabled
        integer, intent(in) :: k_wat_sfc, k_wat_btm, k_sed_sfc
        integer :: i

        if (.not. allocated(self%sources)) return

        do i = 1, size(self%sources)
            if (.not. self%sources(i)%active) cycle

            select case (trim(self%sources(i)%domain))

                case ('interior')
                    ! Interior variables can receive fluxes at water_surface,
                    ! water_bottom, or sediment_surface.

                case ('surface')
                    if (trim(self%sources(i)%target_domain) /= 'water_surface') then
                        error stop 'Surface state variable source "'//trim(self%sources(i)%key)// &
                            '" must use target_domain water_surface.'
                    end if

                case ('bottom')
                    if (sediments_enabled) then
                        error stop 'Bottom state variable source "'//trim(self%sources(i)%key)// &
                            '" cannot be used when sediments are enabled.'
                    end if

                    if (trim(self%sources(i)%target_domain) /= 'water_bottom') then
                        error stop 'Bottom state variable source "'//trim(self%sources(i)%key)// &
                            '" must use target_domain water_bottom.'
                    end if

                case default
                    error stop 'Source "'//trim(self%sources(i)%key)// &
                        '" has unsupported FABM state domain "'//trim(self%sources(i)%domain)//'".'
            end select

            select case (trim(self%sources(i)%target_domain))

                case ('water_surface')
                    self%sources(i)%k_target = k_wat_sfc

                case ('water_bottom')
                    self%sources(i)%k_target = k_wat_btm

                case ('sediment_surface')
                    if (.not. sediments_enabled) then
                        error stop 'Source target_domain sediment_surface requested, but sediments are disabled.'
                    end if
                    self%sources(i)%k_target = k_sed_sfc

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

end module bio_inputs

