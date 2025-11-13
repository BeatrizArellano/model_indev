module grids
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo

  implicit none
  private
  public :: VerticalGrid, build_water_grid, write_vertical_grid, reorder_grid_bottom_first

  type, public :: VerticalGrid
     integer               :: nz = 0             ! number of layers (top=1 … bottom=nz)
     real(rk)              :: depth  = 0.0_rk    ! local water depth [m], positive
     real(rk), allocatable :: z(:)               ! layer centres [m], positive
     real(rk), allocatable :: dz(:)              ! layer thickness [m], positive
     real(rk), allocatable :: z_w(:)             ! layer interfaces [m], positive, size nz+1
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
        type(VerticalGrid),    intent(out) :: grid
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
            grid%depth  = depth
            grid%dz = dz_u
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
                nz = nfull + 1
                call allocate_grid(grid, nz)
                grid%depth = depth
                do i = 1, nfull
                    grid%dz(i) = dz_cfg
                end do
                grid%dz(nz) = rem
                call calculate_layer_depths(grid)
            else
                ! Merge small remainder into the last full cell
                nz = nfull
                call allocate_grid(grid, nz)
                grid%depth = depth
                do i = 1, nz-1
                    grid%dz(i) = dz_cfg
                end do
                grid%dz(nz) = dz_cfg + rem
                call calculate_layer_depths(grid)
            end if

        case default
            call zero_grid(grid)
            if (present(ierr)) then; ierr = 99; return
            else; stop 'build_water_grid: unknown method'
            end if
            end select
    end subroutine build_water_grid
    
    ! Some subroutines require the grids to be inverted
    subroutine reorder_grid_bottom_first(g)
        type(VerticalGrid), intent(inout) :: g
        integer :: N, i, lbz, ubz, lbdz, ubdz
        real(rk), allocatable :: z_new(:), dz_new(:), zw_new(:)
        real(rk) :: depth_check

        N = g%nz
        if (N <= 0) return

        ! ---- Flip centers (1..N) ----
        lbz = lbound(g%z,1);  ubz = ubound(g%z,1)    ! Lowest and upper indices found in the array
        if (ubz-lbz+1 /= N) stop 'reorder_grid_bottom_first: size(z) mismatch'
        allocate(z_new(1:N))
        z_new  = g%z (ubz:lbz:-1)
        !do i = 1, N 
        !    z_new(i) = g%z(lbz + (N - i)) 
        !end do
        
        ! ---- Flip layer thicknesses (1..N) ----
        lbdz = lbound(g%dz,1);  ubdz = ubound(g%dz,1)
        allocate(dz_new(1:N))
        dz_new = g%dz(ubdz:lbdz:-1)
        !do i = 1, N 
        !    dz_new(i) = g%dz(lbdz + (N - i)) 
        !    if (dz_new(i) <= 0.0_rk) stop 'reorder_grid_bottom_first: dz <= 0 after flip' 
        !end do

        ! ---- rebuild interfaces (0..N) from dz (depth-from-surface) ----
        allocate(zw_new(0:N))
        depth_check = sum(dz_new)
        zw_new(0) = depth_check                ! seabed (depth value)
        do i = 1, N
            zw_new(i) = zw_new(i-1) - dz_new(i)  ! steps down to ~0 at surface
        end do

        ! ---- commit ----
        g%z   = z_new
        g%dz  = dz_new
        g%z_w = zw_new
        g%depth = depth_check                  ! keep depth consistent
    end subroutine reorder_grid_bottom_first


    subroutine write_vertical_grid(grid, filename, ierr)
        type(VerticalGrid), intent(in)  :: grid
        character(*),    intent(in)  :: filename
        integer, optional, intent(out) :: ierr

        integer :: u, i
        if (present(ierr)) ierr = 0

        if (grid%nz <= 0 .or. .not. allocated(grid%z) .or. .not. allocated(grid%z_w)) then
            if (present(ierr)) then; ierr = 1; return
            else; stop 'write_vertical_grid: grid not initialized'
            end if
        end if

        ! Get a free unit and open file (overwrite if exists)
        inquire (iolength=u) u
        open(newunit=u, file=trim(filename), status='replace', action='write', form='formatted')

        write(u,'(A)') '# layer  depth_center_m  interface_top_m  interface_bottom_m'
        do i = 1, grid%nz
            write(u,'(I6,1X,F12.6,1X,F12.6,1X,F12.6)') i, grid%z(i), grid%z_w(i), grid%z_w(i+1)
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
        allocate(g%z_w(nz+1))
    end subroutine allocate_grid

    subroutine zero_grid(g)
        type(VerticalGrid), intent(inout) :: g
        if (allocated(g%z))   deallocate(g%z)
        if (allocated(g%dz))  deallocate(g%dz)
        if (allocated(g%z_w)) deallocate(g%z_w)
        g%nz = 0
        g%depth  = 0._rk
    end subroutine zero_grid

    ! Depth convention:
    !   z_w(1)=0 (surface), z_w increases downward to  the bottom
    !   dz(i)  = z_w(i+1) - z_w(i) > 0
    !   z(i)   = 0.5*(z_w(i) + z_w(i+1)) > 0
    subroutine calculate_layer_depths(g)
        type(VerticalGrid), intent(inout) :: g
        integer :: i
        real(rk) :: total

        g%z_w(1) = 0._rk
        do i = 1, g%nz
            g%z_w(i+1) = g%z_w(i) + g%dz(i)
            g%z(i)     = 0.5_rk * (g%z_w(i) + g%z_w(i+1))
        end do

        ! Match the exact depth in water column
        total = g%z_w(g%nz+1)
        if (abs(total - g%depth) > max(1.0e-12_rk, 1.0e-9_rk * g%depth)) then
            g%z_w(g%nz+1) = g%depth
            g%dz(g%nz)    = g%z_w(g%nz+1) - g%z_w(g%nz)
            g%z(g%nz)     = 0.5_rk * (g%z_w(g%nz) + g%z_w(g%nz+1))
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