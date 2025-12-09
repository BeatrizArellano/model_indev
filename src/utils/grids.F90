module grids
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams

  implicit none
  private
  public :: VerticalGrid, build_water_grid, build_sediment_grid, write_vertical_grid

  type, public :: VerticalGrid
    integer               :: nz = 0             ! number of layers (bottom=1 … surface=nz)
    real(rk)              :: depth  = 0.0_rk    ! local water depth [m], positive
    real(rk), allocatable :: z(:)               ! layer centres [m], positive (1=bottom … nz=top)
    real(rk), allocatable :: dz(:)              ! layer thickness [m], positive
    real(rk), allocatable :: z_w(:)             ! layer interfaces [m], 0=bottom … nz=surface
  end type

  ! Defaults for the water column
  real(rk), parameter :: default_dz = 2.0_rk
  real(rk), parameter :: min_dz     = 0.1_rk
  integer,  parameter :: min_nz     = 3

  ! Sediment-grid defaults
  real(rk), parameter :: sed_growth_factor = 2.0_rk     ! geometric growth
  real(rk), parameter :: sed_min_dz        = 1.0e-4_rk  ! minimum allowed thickness for any sediment layer [m]
  real(rk), parameter :: sed_min_depth     = 3.0e-2_rk  ! Minimum depth for sediments [m] (3 cm)
  real(rk), parameter :: sed_default_dzmin = 1.0e-3_rk  ! default top layer [m]  (1 mm)
  real(rk), parameter :: sed_default_dzmax = 2.0e-2_rk  ! default max layer [m]  (2 cm)

contains

    !------------------------------
    ! Build water grid
    !------------------------------
    subroutine build_water_grid(cfg, depth, grid, ierr)
        type(ConfigParams), intent(in)  :: cfg
        real(rk),           intent(in)  :: depth         ! local depth [m], positive
        type(VerticalGrid), intent(out) :: grid
        integer,  optional, intent(out) :: ierr      ! 0 OK; >0 is an error

        character(len=16) :: method
        integer           :: nz_cfg
        real(rk)          :: dz_cfg
        integer           :: nfull, i, nz
        real(rk)          :: rem, dz_u

        if (present(ierr)) ierr = 0

        ! Read configuration choices
        call read_wgrid_config(cfg, nz_cfg, dz_cfg, method)

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
    !   - depth  is total sediment depth [m] below the sediment–water interface.
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
    subroutine build_sediment_grid(cfg, grid, ierr)
        type(ConfigParams), intent(in)  :: cfg
        type(VerticalGrid), intent(out) :: grid
        integer,  optional, intent(out) :: ierr          ! 0 is OK. >0 error codes

        real(rk) :: depth, dz_min, dz_max
        real(rk) :: rem, dz_next, dz_cur
        real(rk), allocatable :: dz_tmp(:)
        real(rk) :: orig_bottom, depth_extra
        real(rk) :: dz
        integer :: nmax, n, i
        logical :: uniform_layers

        if (present(ierr)) ierr = 0

        ! Read user-defined choices for sediment grid
        call read_sgrid_config(cfg, depth, dz_min, dz_max)

        !------------------------------------------------------------------
        ! Special case: uniform grid (dz_min == dz_max)
        !------------------------------------------------------------------
        uniform_layers = (abs(dz_max - dz_min) <= 1.0e-6_rk * dz_min)

        if (uniform_layers) then
            dz = dz_min

            ! Number of layers: increase depth to fit full layers
            n = int(ceiling(depth / dz))
            if (n < 2) n = 2

            call allocate_grid(grid, n)

            ! Bottom-first indexing, thickness is uniform
            grid%dz(:) = dz

            ! Actual sediment depth is now exactly n * dz (>= initial depth)
            grid%depth = dz * n

            call calculate_layer_depths(grid)
            return
        end if

        !------------------------------------------------------------------
        ! General case: grid with layer thicknesses increasing geometrically
        !------------------------------------------------------------------

        ! Maximum number of layers (temporary array): worst case all have a thickness dz_min
        nmax = int(ceiling(depth / max(dz_min, sed_min_dz))) + 10
        allocate(dz_tmp(nmax))

        rem     = depth     ! Remainder depth
        dz_next = dz_min    ! dz for the next layer
        n       = 0         ! Layer number

        ! Build thicknesses from the sediment-water interface downward
        do while (rem > sed_min_dz .and. n < nmax)
            dz_cur = min(dz_next, dz_max)   ! Thickness for the current layer
            n      = n + 1
            ! If what's left after this layer is very small, just make this the last layer
            if (rem - dz_cur <= sed_min_dz) then
                dz_tmp(n) = rem
                rem     = 0._rk
                exit
            else
                dz_tmp(n) = dz_cur
                rem     = rem - dz_cur
                ! Grow geometrically until dz_max is reached
                if (dz_next < dz_max) dz_next = dz_next * sed_growth_factor
            end if
        end do

        ! Ensure at least 2 layers
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

        ! Reverse: dz_tmp(1)=top to grid%dz(n)=top
        do i = 1, n
            grid%dz(i) = dz_tmp(n - i + 1)
        end do

        deallocate(dz_tmp)

        ! --- Enforce a full bottom layer thickness of dz_max ---
        orig_bottom = grid%dz(1)

        if (orig_bottom < dz_max) then
            depth_extra  = dz_max - orig_bottom
            grid%dz(1)   = dz_max
            grid%depth   = depth + depth_extra
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


    subroutine read_wgrid_config(cfg, nz, dz, method)
        type(ConfigParams), intent(in)  :: cfg
        integer,            intent(out) :: nz
        real(rk),           intent(out) :: dz
        character(len=*),   intent(out) :: method   ! 'dz' or 'nlayers'

        logical :: has_dz, has_nz

        ! Safe defaults
        nz     = -1
        dz     = default_dz
        method = 'dz'

        has_dz = .not. cfg%is_disabled('grid.water.dz')
        has_nz = .not. cfg%is_disabled('grid.water.nlayers')

        if (has_nz .and. has_dz) then
            ! Prefer dz when both are provided
            dz     = cfg%get_param_num('grid.water.dz', required=.true., positive=.true., &
                                        min=min_dz, finite=.true.)
            method = 'dz'

        else if (has_nz) then
            nz     = cfg%get_param_int('grid.water.nlayers', min=min_nz)
            method = 'nlayers'

        else if (has_dz) then
            dz     = cfg%get_param_num('grid.water.dz', required=.true., positive=.true., &
                                        min=min_dz, finite=.true.)
            method = 'dz'
        end if
    end subroutine read_wgrid_config

    ! Read parameters from the sediment grid configuration
    subroutine read_sgrid_config(cfg, depth, dz_min, dz_max)
        type(ConfigParams), intent(in)  :: cfg
        real(rk),           intent(out) :: depth   ! sediment depth [m]
        real(rk),           intent(out) :: dz_min  ! min layer thickness at top [m]
        real(rk),           intent(out) :: dz_max  ! max layer thickness [m]


        depth = cfg%get_param_num('grid.sediments.depth', required=.true., positive=.true., &
                                    min=sed_min_depth, finite=.true.)

        dz_min = cfg%get_param_num('grid.sediments.dz_min', positive=.true., &
                                    min=sed_min_dz, default=sed_default_dzmin, finite=.true.)
        
        dz_max = cfg%get_param_num('grid.sediments.dz_max', positive=.true., &
                                    min=sed_min_dz, default=sed_default_dzmax, finite=.true.)

        if (dz_max < dz_min) then
            stop 'read_sgrid_config: dz_max must be >= dz_min in the sediment grid.'
        end if 
        if (dz_min >= depth) then
            stop 'read_sgrid_config: dz_min must be < depth in the sediment grid.'
        end if

        if (dz_max >= depth) then
            stop 'read_sgrid_config: dz_max must be < depth in the sediment grid.'
        end if

    end subroutine read_sgrid_config

end module grids