module trawling

    use bio_types,        only: BioEnv
    use event_types,      only: Event, EventData
    use event_utils,      only: fatal
    use precision_types,  only: rk
    use read_config_yaml, only: ConfigParams
    

    implicit none
    private

    public :: TrawlingData
    public :: trawling_prepare
    public :: trawling_apply
    public :: trawling_clear

    !==================================================================
    ! TrawlingData-specific data
    !==================================================================
    type, extends(EventData) :: TrawlingData
        real(rk) :: erosion_depth       = 0._rk   ! [cm]
        real(rk) :: penetration_depth   = 0._rk   ! total depth affected by the trawling gear [cm]
        real(rk) :: resuspension_height = 0._rk   ! [m]
        ! Indices of relevant layers
        integer :: erosion_idx      = 0           ! index in BE%sed_grid
        integer :: penetration_idx  = 0           ! index in BE%sed_grid
        integer :: resuspension_idx = 0           ! index in BE%wat_grid

        ! Fauna disturbance
        logical  :: affect_fauna    = .false.
        real(rk) :: fauna_depletion = 0._rk       ! Fraction of faunal activity removed [-]


        real(rk) :: erosion_depth_eff
        real(rk) :: penetration_depth_eff
    end type TrawlingData


contains

    !==================================================================
    ! Prepare trawling events
    !==================================================================
    subroutine trawling_prepare(evt, cfg, BE)
        type(Event),        intent(inout) :: evt
        type(ConfigParams), intent(in)    :: cfg
        type(BioEnv),       intent(in)    :: BE

        character(len=:), allocatable :: basekey
        real(rk) :: sediment_depth_cm

        if (allocated(evt%details)) deallocate(evt%details)
        allocate(TrawlingData :: evt%details)

        select type (TR => evt%details)
        type is (TrawlingData)

            basekey = 'events.' // trim(evt%name)

            !----------------------------------------------------------
            ! TrawlingData is an instantaneous event
            !----------------------------------------------------------
            evt%is_instantaneous       = .true.  
            evt%changes_tendency       = .false.            
            evt%changes_state_directly = .false.

            !----------------------------------------------------------
            ! Check requirements
            !----------------------------------------------------------
            ! TrawlingData requires an active sediment domain
            if (.not. BE%params%sediments_enabled) then
                call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                           '" requires sediments to be enabled.')
            end if

            if (BE%nsed <= 0) then
                call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                           '" requires at least one sediment layer.')
            end if

            !----------------------------------------------------------
            ! Read parameters
            !----------------------------------------------------------
            TR%erosion_depth = cfg%get_param_num(basekey//'.erosion_depth', required=.true., finite=.true., min=0._rk)
            TR%penetration_depth = cfg%get_param_num(basekey//'.penetration_depth', required=.true., finite=.true., positive=.true.)
            TR%resuspension_height = cfg%get_param_num(basekey//'.resuspension_height', required=.true., finite=.true., min=0._rk)

            TR%affect_fauna = cfg%get_param_logical(basekey//'.affect_fauna', default=.false.)
            if (TR%affect_fauna) then
                TR%fauna_depletion = cfg%get_param_num(basekey//'.fauna_depletion', required=.true., finite=.true., min=0._rk, max=1._rk)
            else
                TR%fauna_depletion = 0._rk
            end if


            !----------------------------------------------------------
            ! Validations
            !----------------------------------------------------------
            sediment_depth_cm = 100._rk * BE%sed_grid%depth

            if (TR%erosion_depth > TR%penetration_depth) then
                call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                           '": erosion_depth cannot exceed the penetration depth of trawling device.')
            end if

            if (TR%penetration_depth > sediment_depth_cm) then
                call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                           '": penetration_depth exceeds the sediment depth.')
            end if

            if (TR%affect_fauna .and. .not. BE%SED%use_bioturbation .and. .not. BE%SED%use_bioirrigation) then
                call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                           '": affect_fauna is enabled, but both bioturbation and bioirrigation are off.')
            end if

            if (TR%affect_fauna) then

                if (BE%SED%use_bioturbation .and. .not. BE%SED%output_bioturb_dynamic) then
                    call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                               '": affect_fauna requires bioturbation_mode = dynamic.')
                end if

                if (BE%SED%use_bioirrigation .and. .not. BE%SED%output_bioirr_dynamic) then
                    call fatal('trawling:prepare', 'TrawlingData event "'//trim(evt%name)// &
                               '": affect_fauna requires bioirrigation_mode = dynamic.')
                end if

            end if

            !---------------------------------------
            ! Find relevant indices
            !---------------------------------------
            TR%erosion_idx = find_deepest_full_layer(BE%sed_grid%z_w, 0.01_rk * TR%erosion_depth, TR%erosion_depth_eff)
            TR%penetration_idx = find_deepest_full_layer(BE%sed_grid%z_w, 0.01_rk * TR%penetration_depth, TR%penetration_depth_eff)

            if (TR%penetration_idx == 0) then
                ! Minimum penetration depth: top sediment layer
                TR%penetration_idx = BE%nsed
                TR%penetration_depth_eff = BE%sed_grid%dz(BE%nsed)
            end if

            if (TR%resuspension_height <= 0._rk) then
                ! Bottom water layer only
                TR%resuspension_idx = 1

            else if (TR%resuspension_height >= BE%wat_grid%depth) then
                ! Entire water column
                TR%resuspension_idx = BE%wat_grid%nz
            else
                TR%resuspension_idx = find_shallowest_full_layer_from_bottom(BE%wat_grid%z_w, TR%resuspension_height)   
                if (TR%resuspension_idx == 0) TR%resuspension_idx = 1
            end if

        end select

    end subroutine trawling_prepare


    !==================================================================
    ! Apply instantaneous trawling disturbance
    !==================================================================
    subroutine trawling_apply(evt, BE)

        type(Event),  intent(in)    :: evt
        type(BioEnv), intent(inout) :: BE

        integer  :: ivar
        integer  :: wat_shallow_full_idx

        real(rk) :: mass_eroded
        real(rk) :: mass_penetration
        real(rk) :: mass_remaining

        real(rk) :: solid_capacity, porewater_capacity
        real(rk) :: mass_in_sediment, mass_in_water
        real(rk) :: water_thickness
        real(rk) :: mixed_conc
        real(rk) :: delta_conc

        ! Debugging variables to check mass conservation
        !real(rk), allocatable :: mass_before(:), mass_after(:)

        select type (TR => evt%details)
        type is (TrawlingData)

            !---------Mass conservation check (Debug)------
            !allocate(mass_before(BE%BS%n_interior))
            !allocate(mass_after(BE%BS%n_interior))
            !
            !do ivar = 1, BE%BS%n_interior
            !    mass_before(ivar) = total_tracer_inventory(BE, ivar)
            !end do
            !------------------------------------------------------

            !----------------------------------------------------------
            ! Grid geometry
            !----------------------------------------------------------
            !
            ! Sediment indices are already full-column indices: 1:BE%nsed
            ! Water indices stored in TR are local to BE%wat_grid.
            ! Convert the shallowest affected water index to the
            ! corresponding full-column index.
            !
            wat_shallow_full_idx = BE%k_wat_btm + TR%resuspension_idx - 1

            ! Total water thickness receiving resuspended material
            water_thickness = sum(BE%wat_grid%dz(1:TR%resuspension_idx))

            if (water_thickness <= 0._rk) then
                call fatal('trawling:apply', &
                    'Non-positive water thickness for resuspension.')
            end if

            !----------------------------------------------------------
            ! Particulate tracers
            !----------------------------------------------------------
            do ivar = 1, BE%BS%n_interior

                if (.not. BE%tracer_info(ivar)%is_particulate) cycle

                !------------------------------------------------------
                ! 1. Total particulate inventory within the
                !    penetration layer before disturbance.
                !
                ! Concentrations in sediments are expressed per
                ! solid volume:
                !     inventory = C * (1-phi) * dz
                !
                ! solid_thickness already contains (1-phi)*dz.
                !------------------------------------------------------
                mass_penetration = sum(BE%BS%interior_state(TR%penetration_idx:BE%nsed, ivar) * &
                                       BE%SED%solid_thickness(TR%penetration_idx:BE%nsed))

                !------------------------------------------------------
                ! 2. Inventory resuspended from the eroded layer.
                !
                ! erosion_idx = 0 means that the requested erosion
                ! depth did not contain a complete sediment layer.
                !------------------------------------------------------
                mass_eroded = 0._rk

                if (TR%erosion_idx > 0) then
                    mass_eroded = sum(BE%BS%interior_state(TR%erosion_idx:BE%nsed, ivar) * &
                                      BE%SED%solid_thickness(TR%erosion_idx:BE%nsed))
                end if

                !------------------------------------------------------
                ! 3. Resuspend particulate matter
                !
                ! Water concentrations are expressed per water volume,
                ! so dividing areal inventory by water thickness gives
                ! the concentration increment.
                !------------------------------------------------------
                if (TR%erosion_idx > 0) then
                    delta_conc = mass_eroded / water_thickness
                    BE%BS%interior_state(BE%k_wat_btm:wat_shallow_full_idx, ivar) = &
                                        BE%BS%interior_state(BE%k_wat_btm:wat_shallow_full_idx, ivar) + delta_conc
                end if

                !------------------------------------------------------
                ! 4. Particulate inventory remaining within the
                !    penetration layer after resuspension.
                !------------------------------------------------------
                mass_remaining = mass_penetration - mass_eroded

                !------------------------------------------------------
                ! 5. Homogenising the remaining particulate inventory
                !    throughout the penetration layer.
                !
                ! Because concentrations are per solid volume, divide
                ! by the total solid storage capacity rather than by
                ! the geometrical sediment thickness.
                !------------------------------------------------------
                solid_capacity = sum(BE%SED%solid_thickness(TR%penetration_idx:BE%nsed))
                mixed_conc = mass_remaining / solid_capacity
                BE%BS%interior_state(TR%penetration_idx:BE%nsed, ivar) = mixed_conc        

            end do

            !----------------------------------------------------------
            ! Solute tracers
            !----------------------------------------------------------
            do ivar = 1, BE%BS%n_interior

                if (.not. BE%tracer_info(ivar)%is_solute) cycle

                !------------------------------------------------------
                ! 1. Computing solute inventory in disturbed porewater
                !
                ! Sediment solute concentrations are expressed per
                ! porewater volume:
                !     inventory = C * phi * dz
                ! porewat_thickness already contains phi*dz.
                !------------------------------------------------------
                mass_in_sediment = sum(BE%BS%interior_state(TR%penetration_idx:BE%nsed, ivar) * &
                                    BE%SED%porewat_thickness(TR%penetration_idx:BE%nsed))

                !------------------------------------------------------
                ! 2. Solute inventory in bottom water layer
                !------------------------------------------------------
                mass_in_water = BE%BS%interior_state(BE%k_wat_btm, ivar) * BE%wat_grid%dz(1)

                !------------------------------------------------------
                ! 3. Total porewater volume
                !------------------------------------------------------
                porewater_capacity = sum(BE%SED%porewat_thickness(TR%penetration_idx:BE%nsed))

                !------------------------------------------------------
                ! 4. Homogeneous concentration across porewater & overlying water
                !------------------------------------------------------
                mixed_conc = (mass_in_sediment + mass_in_water) / (porewater_capacity + BE%wat_grid%dz(1))

                !------------------------------------------------------
                ! 5. Apply the same concentration throughout the
                !    disturbed sediment porewater
                !------------------------------------------------------
                BE%BS%interior_state(TR%penetration_idx:BE%nsed, ivar) = mixed_conc

                !------------------------------------------------------
                ! 6. Apply the same concentration in bottom water
                !------------------------------------------------------
                BE%BS%interior_state(BE%k_wat_btm, ivar) = mixed_conc

            end do

            !----------------------------------------------------------
            ! Apply disturbance to benthic-faunal activity
            !----------------------------------------------------------
            if (TR%affect_fauna .and. TR%fauna_depletion > 0.0_rk) then
                BE%SED%faunal_activity = max(0.0_rk, BE%SED%faunal_activity * (1.0_rk - TR%fauna_depletion))
            end if

            !---------Mass conservation check (Debug)------
            !do ivar = 1, BE%BS%n_interior
            !    mass_after(ivar) = total_tracer_inventory(BE, ivar)
            !end do

            !call check_mass_conservation(BE, mass_before, mass_after)
            !deallocate(mass_before, mass_after)
            !-------------------------------------------------------

        class default
            call fatal('trawling:apply', &
                'Invalid event data for trawling event "'// &
                trim(evt%name)//'".')
        end select

    end subroutine trawling_apply


    !==================================================================
    ! Clear event-specific allocated data
    !==================================================================
    subroutine trawling_clear(TR)

        type(TrawlingData), intent(inout) :: TR

        ! Nothing to clear yet.
        ! Add deallocation here if TrawlingData later owns allocatable data.

    end subroutine trawling_clear

    !=================================
    !   Internal
    !=================================
    ! Returns the deepest layer that is fully included between the SWI and
    ! the requested depth, without exceeding that depth.
    !
    ! Grid convention:
    !   z_w(nz) = 0       ! SWI
    !   z_w(0)  = depth   ! bottom of sediment grid
    !   layers are indexed 1:nz from bottom to surface
    !
    ! Returns idx=0 if the requested depth does not include one complete layer.
    !
    ! Optionally returns the effective depth, snapped to the corresponding
    ! sediment interface.
    integer function find_deepest_full_layer(z_w, depth, effective_depth) result(idx)
        real(rk), intent(in)  :: z_w(0:)
        real(rk), intent(in)  :: depth
        real(rk), intent(out), optional :: effective_depth

        integer  :: interface_idx, nz
        real(rk) :: depth_eff

        nz = ubound(z_w, 1)

        idx       = 0
        depth_eff = 0._rk

        ! Move downward from the SWI through sediment interfaces.
        ! Keep accepting complete layers while the next interface
        ! does not exceed the requested depth.
        do interface_idx = nz - 1, 0, -1

            if (z_w(interface_idx) <= depth) then
                idx       = interface_idx + 1
                depth_eff = z_w(interface_idx)
            else
                exit
            end if

        end do

        if (present(effective_depth)) effective_depth = depth_eff

    end function find_deepest_full_layer

    ! Returns the shallowest water layer fully included within a specified
    ! height above the seabed, without exceeding that height.
    !
    ! Grid convention:
    !   z_w(0)  = water depth   ! seabed
    !   z_w(nz) = 0             ! sea surface
    !   layers are indexed 1:nz from bottom to surface
    !
    ! Returns idx=0 if the requested height is smaller than the bottom
    ! water-layer thickness.
    integer function find_shallowest_full_layer_from_bottom(z_w, height, effective_height) result(idx)
        real(rk), intent(in)  :: z_w(0:)
        real(rk), intent(in)  :: height
        real(rk), intent(out), optional :: effective_height

        integer  :: interface_idx, nz
        real(rk) :: height_eff, current_height

        nz = ubound(z_w, 1)

        idx        = 0
        height_eff = 0._rk

        ! Move upward from the seabed through water-layer interfaces.
        do interface_idx = 1, nz

            current_height = z_w(0) - z_w(interface_idx)

            if (current_height <= height) then
                idx        = interface_idx
                height_eff = current_height
            else
                exit
            end if

        end do

        if (present(effective_height)) effective_height = height_eff

    end function find_shallowest_full_layer_from_bottom


    real(rk) function total_tracer_inventory(BE, ivar) result(total_mass)

        type(BioEnv), intent(in) :: BE
        integer,      intent(in) :: ivar

        total_mass = 0._rk

        ! Water-column inventory
        total_mass = total_mass + sum( &
            BE%BS%interior_state(BE%k_wat_btm:BE%k_wat_sfc, ivar) * &
            BE%wat_grid%dz )

        ! Sediment inventory
        if (BE%tracer_info(ivar)%is_particulate) then

            total_mass = total_mass + sum( &
                BE%BS%interior_state(1:BE%nsed, ivar) * &
                BE%SED%solid_thickness )

        else if (BE%tracer_info(ivar)%is_solute) then

            total_mass = total_mass + sum( &
                BE%BS%interior_state(1:BE%nsed, ivar) * &
                BE%SED%porewat_thickness )

        end if

    end function total_tracer_inventory

    subroutine check_mass_conservation(BE, before, after)

        type(BioEnv), intent(in) :: BE
        real(rk),     intent(in) :: before(:), after(:)

        integer  :: ivar
        real(rk) :: abs_error, rel_error
        logical  :: conserved

        conserved = .true.

        do ivar = 1, size(before)

            abs_error = after(ivar) - before(ivar)

            if (abs(before(ivar)) > tiny(1._rk)) then
                rel_error = abs_error / abs(before(ivar))
            else
                rel_error = abs_error
            end if

            if (abs(rel_error) > 1.0e-10_rk) then
                conserved = .false.

                write(*,'(A,A)') &
                    'WARNING: trawling mass conservation error for tracer ', &
                    trim(BE%model%interior_state_variables(ivar)%name)

                write(*,'(A,ES14.6)') '  before = ', before(ivar)
                write(*,'(A,ES14.6)') '  after  = ', after(ivar)
                write(*,'(A,ES14.6)') '  diff   = ', abs_error
                write(*,'(A,ES14.6)') '  rel    = ', rel_error
            end if

        end do

        if (conserved) then
            write(*,*) 'TrawlingData mass conservation check: OK'
        end if

    end subroutine check_mass_conservation

end module trawling