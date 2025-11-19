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


  real(rk), parameter :: default_dz = 2.0
  real(rk), parameter :: min_dz     = 0.2
  integer,  parameter :: min_nz     = 3

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