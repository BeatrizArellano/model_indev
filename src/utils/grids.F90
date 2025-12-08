module grids
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo

  implicit none
  private
  public :: VerticalGrid, build_water_grid, write_vertical_grid

  type, public :: VerticalGrid
    integer               :: nz = 0             ! number of layers (bottom=1 … surface=nz)
    real(rk)              :: depth  = 0.0_rk    ! local water depth [m], positive
    real(rk), allocatable :: z(:)               ! layer centres [m], positive (1=bottom … nz=top)
    real(rk), allocatable :: dz(:)              ! layer thickness [m], positive
    real(rk), allocatable :: z_w(:)             ! layer interfaces [m], 0=bottom … nz=surface
  end type

  ! Defaults for the water column
  real(rk), parameter :: default_dz = 2.0
  real(rk), parameter :: min_dz     = 0.1
  integer,  parameter :: min_nz     = 3

  ! Sediment-grid defaults
  real(rk), parameter :: sed_growth_factor = 2.0_rk    ! geometric growth
  real(rk), parameter :: sed_min_dz       = 1.0e-4_rk  ! minimum allowed thickness [m]
  real(rk), parameter :: sed_default_dzmin = 1.0e-3_rk ! default top layer [m]  (1 mm)
  real(rk), parameter :: sed_default_dzmax = 2.0e-2_rk ! default max layer [m]  (2 cm)

contains

    !------------------------------
    ! Build water grid
    !------------------------------
    subroutine build_water_grid(params, depth, grid, ierr)
        type(ConfigParams), intent(in)  :: params
        real(rk),           intent(in)  :: depth         ! local depth [m], positive
        type(VerticalGrid), intent(out) :: grid
        integer,  optional, intent(out) :: ierr      ! 0 OK; 1 bad depth; 2 nz too thin; 3 dz too thin

        character(len=16) :: method
        integer           :: nz_cfg
        real(rk)          :: dz_cfg
        integer           :: nfull, i, nz
        real(rk)          :: rem, dz_u

        if (present(ierr)) ierr = 0

        ! Read YAML choices
        call read_wgrid_config(params, nz_cfg, dz_cfg, method)

        select case (trim(method))

        case ('nlayers')

            dz_u = depth / real(nz_cfg, rk)
            if (dz_u < min_dz) then
                call zero_grid(grid)
                if (present(ierr)) then; ierr = 2; return
                else; stop 'build_water_grid: depth/nlayers < min_dz'
                end if
            end if

            call allocate_grid(grid, nz_cfg)
            grid%depth = depth
            grid%dz    = dz_u          ! all layers same thickness; bottom→top is irrelevant here
            call calculate_layer_depths(grid)

        case ('dz')

            nfull = int(floor(depth / dz_cfg))
            rem   = depth - real(nfull, rk) * dz_cfg

            ! Ensure at least min_nz layers
            if (nfull < min_nz) then
                call zero_grid(grid)
                if (present(ierr)) then; ierr = 2; return
                else; stop 'build_water_grid: depth/dz < 3'
                end if
            end if

            if (rem >= min_dz) then
                ! One remainder cell + nfull full cells
                nz = nfull + 1
                call allocate_grid(grid, nz)
                grid%depth = depth

                ! Bottom-first convention: put the remainder at the bottom (i=1)
                grid%dz(1) = rem
                do i = 2, nz
                    grid%dz(i) = dz_cfg
                end do

                call calculate_layer_depths(grid)

            else
                ! Merge small remainder into the last full cell
                nz = nfull
                call allocate_grid(grid, nz)
                grid%depth = depth

                ! Bottom-first convention: make the bottom cell thicker
                grid%dz(1) = dz_cfg + rem
                do i = 2, nz
                    grid%dz(i) = dz_cfg
                end do

                call calculate_layer_depths(grid)
            end if

        case default
            call zero_grid(grid)
            if (present(ierr)) then
                ierr = 99; return
            else
                stop 'build_water_grid: unknown method'
            end if
        end select
    end subroutine build_water_grid

    !------------------------------
    ! Build sediment grid
    !   - depth > 0 is total sediment depth [m] below the sediment–water interface.
    !   - z_w(0)    = depth       (bottom of sediment column)
    !   - z_w(nz)   = 0           (sediment–water interface)
    !   - dz(i)     = z_w(i-1) - z_w(i) > 0, bottom to top
    !   - z(i)      = 0.5*(z_w(i-1) + z_w(i)), 1=bottom … nz=top (interface)
    !
    ! Layer thicknesses are constructed TOP-DOWN as a geometric sequence:
    !   dz_top = dz_min, then dz_top*growth, dz_top*growth^2, ...
    !   until dz reaches dz_max; below that all layers have dz = dz_max.
    !
    ! Then the sequence is reversed to match bottom-first indexing.
    ! The bottom cell is finally adjusted to dz_max (if smaller),
    ! and grid%depth is increased accordingly so that:
    !   sum(dz) = grid%depth
    !------------------------------
    subroutine build_sediment_grid(params, depth, grid, ierr)
        type(ConfigParams), intent(in)  :: params
        real(rk),           intent(in)  :: depth         ! requested sediment depth [m], positive
        type(VerticalGrid), intent(out) :: grid
        integer,  optional, intent(out) :: ierr          ! 0 OK; >0 error codes

        real(rk) :: dz_min_sed, dz_max_sed
        real(rk) :: rem, dz_cur, dz_eff, tol
        real(rk), allocatable :: dz_tmp(:)
        real(rk) :: orig_bottom, depth_extra
        integer :: nmax, n, i

        if (present(ierr)) ierr = 0

        if (depth <= 0._rk) then
            call zero_grid(grid)
            if (present(ierr)) then
                ierr = 1
                return
            else
                stop 'build_sediment_grid: depth must be positive'
            end if
        end if

        ! Read YAML choices for sediment grid
        call read_sgrid_config(params, dz_min_sed, dz_max_sed)

        ! Safety checks
        if (dz_min_sed <= 0._rk .or. dz_max_sed <= 0._rk) then
            call zero_grid(grid)
            if (present(ierr)) then
                ierr = 2
                return
            else
                stop 'build_sediment_grid: dz_min_sed and dz_max_sed must be positive'
            end if
        end if
        if (dz_min_sed > dz_max_sed) then
            ! Swap if user inverted them
            tol        = dz_min_sed
            dz_min_sed = dz_max_sed
            dz_max_sed = tol
        end if

        ! Upper bound on number of layers: worst case all at dz_min_sed
        nmax = int(ceiling(depth / max(dz_min_sed, sed_min_dz))) + 10
        if (nmax < 2) nmax = 2
        allocate(dz_tmp(nmax))

        rem    = depth
        dz_cur = dz_min_sed
        n      = 0

        ! Build thicknesses from TOP (interface) downward
        do while (rem > sed_min_dz .and. n < nmax)
            dz_eff = min(dz_cur, dz_max_sed)

            ! If what's left after this cell would be very small, just make this the last cell
            if (rem - dz_eff <= sed_min_dz) then
                n       = n + 1
                dz_tmp(n) = rem
                rem     = 0._rk
                exit
            else
                n       = n + 1
                dz_tmp(n) = dz_eff
                rem     = rem - dz_eff

                ! Grow geometrically until dz_max_sed is reached
                if (dz_cur < dz_max_sed) dz_cur = dz_cur * sed_growth_factor
            end if
        end do

        ! If we ran out of loop but still have significant remainder, append it
        if (rem > sed_min_dz) then
            if (n == nmax) then
                call zero_grid(grid)
                if (present(ierr)) then
                    ierr = 3
                    return
                else
                    stop 'build_sediment_grid: exceeded maximum number of sediment layers'
                end if
            end if
            n       = n + 1
            dz_tmp(n) = rem
            rem     = 0._rk
        end if

        ! Ensure at least 2 layers for a meaningful vertical profile
        if (n < 2) then
            ! Just split the depth into two equal layers
            n         = 2
            dz_tmp(1) = 0.5_rk * depth
            dz_tmp(2) = 0.5_rk * depth
        end if

        ! Now allocate the final grid with n layers (bottom-first indexing)
        call allocate_grid(grid, n)

        ! Requested depth so far (may be increased below)
        grid%depth = depth

        ! Reverse: dz_tmp(1)=top → grid%dz(n)=top
        do i = 1, n
            grid%dz(i) = dz_tmp(n - i + 1)
        end do

        deallocate(dz_tmp)

        ! --- Enforce a full bottom layer thickness of dz_max_sed if needed ---
        orig_bottom = grid%dz(1)

        if (orig_bottom < dz_max_sed) then
            depth_extra  = dz_max_sed - orig_bottom
            grid%dz(1)   = dz_max_sed
            grid%depth   = depth + depth_extra
        else
            ! No change: keep original bottom thickness and requested depth
            grid%depth = depth
        end if

        ! Compute z, z_w using existing helper (bottom-first)
        call calculate_layer_depths(grid)
    end subroutine build_sediment_grid


    subroutine write_vertical_grid(grid, filename, ierr)
        type(VerticalGrid), intent(in)  :: grid
        character(*),       intent(in)  :: filename
        integer,  optional, intent(out) :: ierr

        integer :: u, i
        if (present(ierr)) ierr = 0

        if (grid%nz <= 0 .or. .not. allocated(grid%z) .or. .not. allocated(grid%z_w)) then
            if (present(ierr)) then
                ierr = 1; return
            else
                stop 'write_vertical_grid: grid not initialized'
            end if
        end if

        inquire (iolength=u) u
        open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

        write(u,'(A)') 'layer  depth_center_m  thickness_m  interface_top_m  interface_bottom_m'
        do i = grid%nz, 1,-1
            write(u,'(I6,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6)') &
                i, grid%z(i), grid%dz(i), grid%z_w(i), grid%z_w(i-1)
        end do

        close(u)
    end subroutine write_vertical_grid



    !------------------------------
    ! Helpers
    !------------------------------

    subroutine allocate_grid(g, nz)
        type(VerticalGrid), intent(inout) :: g
        integer,            intent(in)    :: nz
        g%nz = nz
        if (allocated(g%z))   deallocate(g%z)
        if (allocated(g%dz))  deallocate(g%dz)
        if (allocated(g%z_w)) deallocate(g%z_w)
        allocate(g%z(nz))
        allocate(g%dz(nz))
        allocate(g%z_w(0:nz))     ! 0..nz: 0=bottom, nz=surface
    end subroutine allocate_grid


    subroutine zero_grid(g)
        type(VerticalGrid), intent(inout) :: g
        if (allocated(g%z))   deallocate(g%z)
        if (allocated(g%dz))  deallocate(g%dz)
        if (allocated(g%z_w)) deallocate(g%z_w)
        g%nz = 0
        g%depth  = 0._rk
    end subroutine zero_grid

    ! Depth convention (bottom-first indexing):
    !   z_w(0)    = bottom interface depth (= depth > 0)
    !   z_w(nz)   = surface (0)
    !   dz(i)     = z_w(i-1) - z_w(i) > 0, layer i thickness, bottom→top
    !   z(i)      = 0.5*(z_w(i-1) + z_w(i)), layer centres, 1=bottom … nz=top
    subroutine calculate_layer_depths(g)
        type(VerticalGrid), intent(inout) :: g
        integer :: i
        real(rk) :: tol

        if (g%nz <= 0) return

        ! Set bottom interface exactly at depth
        g%z_w(0) = g%depth

        ! Depth for layer interfaces and centres from the bottom to the surface
        do i = 1, g%nz
            g%z_w(i) = g%z_w(i-1) - g%dz(i)              ! move up (shallower)
            g%z(i)   = 0.5_rk * (g%z_w(i-1) + g%z_w(i))  ! centre depth
        end do

        ! Tolerance for matching the surface at 0 m
        tol = max(1.0e-12_rk, 1.0e-9_rk * g%depth)

        ! Adjust the top cell so the surface is exactly at 0
        if (abs(g%z_w(g%nz)) > tol) then
            g%z_w(g%nz) = 0._rk
            g%dz(g%nz)  = g%z_w(g%nz-1) - g%z_w(g%nz)
            g%z(g%nz)   = 0.5_rk * (g%z_w(g%nz-1) + g%z_w(g%nz))
        end if
    end subroutine calculate_layer_depths


    subroutine read_wgrid_config(params, nz, dz, method)
        type(ConfigParams), intent(in)  :: params
        integer,            intent(out) :: nz
        real(rk),           intent(out) :: dz
        character(len=*),   intent(out) :: method   ! 'dz' or 'nlayers'

        logical :: has_dz, has_nz

        ! Safe defaults
        nz     = -1
        dz     = default_dz
        method = 'dz'

        has_dz = .not. params%is_disabled('grid.water.dz')
        has_nz = .not. params%is_disabled('grid.water.nlayers')

        if (has_nz .and. has_dz) then
            ! Prefer dz when both are provided
            dz     = params%get_param_num('grid.water.dz', required=.true., positive=.true., &
                                        min=min_dz, finite=.true.)
            method = 'dz'

        else if (has_nz) then
            nz     = params%get_param_int('grid.water.nlayers', min=min_nz)
            method = 'nlayers'

        else if (has_dz) then
            dz     = params%get_param_num('grid.water.dz', required=.true., positive=.true., &
                                        min=min_dz, finite=.true.)
            method = 'dz'
        end if
    end subroutine read_wgrid_config

end module grids