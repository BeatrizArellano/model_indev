! src/io/output_writer.F90
module output_writer
  use precision_types,  only: rk, lk
  use netcdf,           only: NF90_UNLIMITED
  use netcdf_io,        only: NcFile, nc_create, nc_redef, nc_enddef, nc_close, nc_sync, &
                              nc_def_dim, nc_def_var_real, nc_def_var_double, &
                              nc_put_att_str, nc_put_att_real, nc_write_real
  use grids,            only: VerticalGrid
  use variable_registry, only: VarMetadata

  implicit none
  private
  public :: OutputWriter


  type :: VarID
    character(len=:), allocatable :: name   ! NetCDF variable name
    integer :: n_space_dims = 0
    character(len=16) :: vert_coord = 'none'  ! 'centre','interface','none',...
    integer :: col_index = 0                  ! which column in the aggregated arrays
  end type VarID


  !========================
  ! Writer 
  !========================
  type :: OutputWriter
     type(NcFile) :: db
     character(:), allocatable :: filename
     character(:), allocatable :: interval_statistic  ! 'mean' | 'instant'
     character(:), allocatable :: time_units          ! "seconds since YYYY-MM-DD hh:mm:ss"
     character(:), allocatable :: calendar_name       ! CF calendar value

     integer(lk) :: N  = 0_lk   ! centers
     integer(lk) :: Ni = 0_lk   ! interfaces (= N+1)
     integer(lk) :: time_index = 0_lk

     ! Counters for variable types
     integer :: n_centre = 0      ! number of center-profile vars
     integer :: n_iface  = 0      ! number of interface-profile vars
     integer :: n_scalar = 0      ! number of scalar vars

     type(VarID), allocatable :: centre_varids(:) ! Ids for variables reported at layers' centres
     type(VarID), allocatable :: iface_varids(:)
     type(VarID), allocatable :: scalar_varids(:)

     ! Working arrays
     real(rk), allocatable :: center_row(:)  ! size N
     real(rk), allocatable :: iface_row(:)   ! size Ni
  contains
     procedure :: open_file    => outputwriter_open_file
     procedure :: append_record=> outputwriter_append_record
     procedure :: close_file   => outputwriter_close_file
     procedure :: sync_file     => outputwriter_sync_file   
  end type OutputWriter

contains  

  !========================
  ! Open file, define coordinates and set attrubutes
  !========================
  subroutine outputwriter_open_file(this, path, grid, interval_statistic, time_units, &
                                    calendar_name, vars, title, n_sed)
    class(OutputWriter), intent(inout) :: this
    character(*),        intent(in)    :: path
    type(VerticalGrid),  intent(in)    :: grid        ! bottom→surface order internally
    character(*),        intent(in)    :: interval_statistic
    character(*),        intent(in)    :: time_units       ! "seconds since ..."
    character(*),        intent(in)    :: calendar_name
    type(VarMetadata),   intent(in)    :: vars(:)
    character(*),        intent(in), optional :: title
    integer,             intent(in), optional :: n_sed

    integer :: dim_t, dim_z, dim_zw
    integer :: j, ic, ii, is 
    integer :: var_time, var_z, var_zw, var_mask, varid_dummy
    integer :: nsed_local
    character(len=:), allocatable :: cm
    real(rk), allocatable :: z_out(:), zw_out(:)
    real(rk), allocatable :: mask_out(:)

    ! sizes
    this%N  = grid%nz
    this%Ni = grid%nz + 1

    nsed_local = 0
    if (present(n_sed)) nsed_local = max(0, n_sed)
    if (nsed_local > int(this%N)) error stop 'outputwriter_open_file: n_sed > N'

    this%filename           = trim(path)
    this%interval_statistic = trim(interval_statistic)
    this%time_units         = trim(time_units)
    this%calendar_name      = trim(calendar_name)
    this%time_index         = 0_lk


    ! Create and define coordinates
    call nc_create(this%db, this%filename, overwrite=.true.)
    !call nc_redef(this%db) <- Fails with old compilers

    dim_t  = nc_def_dim(this%db, 'time',  NF90_UNLIMITED)
    dim_z  = nc_def_dim(this%db, 'depth',     int(this%N, kind=4))
    dim_zw = nc_def_dim(this%db, 'depth_interface',   int(this%Ni, kind=4))

    var_time = nc_def_var_double(this%db, 'time', [dim_t])              ! NOTE: uses rk; prefer rk=real64
    var_z    = nc_def_var_real(this%db, 'depth',    [dim_z])
    var_zw   = nc_def_var_real(this%db, 'depth_interface',  [dim_zw])

    var_mask = nc_def_var_real(this%db, 'sediment_mask', [dim_z])
    call nc_put_att_str(this%db, 'sediment_mask', 'long_name', 'Mask for sediment layers (1=sediment, 0=water)')
    call nc_put_att_str(this%db, 'sediment_mask', 'units', '1')
    call nc_put_att_str(this%db, 'sediment_mask', 'flag_values', '0 1')
    call nc_put_att_str(this%db, 'sediment_mask', 'flag_meanings', 'water sediment')
    call nc_put_att_str(this%db, 'sediment_mask', 'coordinates', 'depth')

    cm = cf_cell_methods_from_stat(interval_statistic) ! Write the statistic following the CF-metadata convention to set as attribute

    this%n_centre = 0
    this%n_iface  = 0
    this%n_scalar = 0

    do j = 1, size(vars)    ! Count variables that are at the centres, interfaces and scalars
        if (.not. vars(j)%output) cycle

        select case (vars(j)%n_space_dims)
          case (1)
              select case (trim(vars(j)%vert_coord))
                case ('centre')
                    this%n_centre = this%n_centre + 1
                case ('interface')
                    this%n_iface  = this%n_iface  + 1
                case default
                    error stop 'outputwriter_open_file: unsupported vert_coord for profile.'
                end select
          case (0)
              this%n_scalar = this%n_scalar + 1
        case default
            error stop 'outputwriter_open_file: unsupported n_space_dims.'
        end select
    end do

    ! Allocate VarID arrays and working rows
    if (this%n_centre > 0) then
        allocate(this%centre_varids(this%n_centre))
        allocate(this%center_row(this%N))
    end if

    if (this%n_iface > 0) then
        allocate(this%iface_varids(this%n_iface))
        allocate(this%iface_row(this%Ni))
    end if

    if (this%n_scalar > 0) then
        allocate(this%scalar_varids(this%n_scalar))
    end if

    ic = 0
    ii = 0
    is = 0

    ! === Organise variables by needed coordinates (at layers' centres, interfaces or scalars)
    do j = 1, size(vars)
        if (.not. vars(j)%output) cycle
        select case (vars(j)%n_space_dims)
          case (1)
              select case (trim(vars(j)%vert_coord))
                case ('centre')                     ! Variables at layers' centres
                    ! time x z (cell centers)                    
                    ic = ic + 1                    
                    this%centre_varids(ic)%name        = trim(vars(j)%name)
                    this%centre_varids(ic)%n_space_dims = vars(j)%n_space_dims
                    this%centre_varids(ic)%vert_coord   = trim(vars(j)%vert_coord)
                    this%centre_varids(ic)%col_index    = ic   ! column in centre_data(:, ic)
                    varid_dummy = nc_def_var_real(this%db, trim(vars(j)%name), [dim_t, dim_z])
                    
                case ('interface')
                    ! time x z_w (interfaces)              
                    ii = ii + 1                    
                    this%iface_varids(ii)%name         = trim(vars(j)%name)
                    this%iface_varids(ii)%n_space_dims = vars(j)%n_space_dims
                    this%iface_varids(ii)%vert_coord   = trim(vars(j)%vert_coord)
                    this%iface_varids(ii)%col_index    = ii    ! column in iface_data(:, ii)
                    varid_dummy = nc_def_var_real(this%db, trim(vars(j)%name), [dim_t, dim_zw])
                case default
                    error stop 'outputwriter_open_file: unsupported vert_coord for profile.'
                end select

          case (0)
              ! scalar: time only
              is = is + 1
              this%scalar_varids(is)%name         = trim(vars(j)%name)
              this%scalar_varids(is)%n_space_dims = vars(j)%n_space_dims
              this%scalar_varids(is)%vert_coord   = trim(vars(j)%vert_coord)
              this%scalar_varids(is)%col_index    = is       ! index in scalar_data(is)
              varid_dummy = nc_def_var_real(this%db, trim(vars(j)%name), [dim_t])

          case default
              error stop 'outputwriter_open_file: unsupported n_space_dims.'
        end select        

        ! Set optional attributes 
        if (allocated(vars(j)%standard_name)) then
            if (len_trim(vars(j)%standard_name) > 0) &
                call nc_put_att_str(this%db, trim(vars(j)%name), 'standard_name', trim(vars(j)%standard_name))
        end if
        if (allocated(vars(j)%long_name)) then
            if (len_trim(vars(j)%long_name) > 0) &
                call nc_put_att_str(this%db, trim(vars(j)%name), 'long_name', trim(vars(j)%long_name))
        end if
        if (allocated(vars(j)%units)) then
            if (len_trim(vars(j)%units) > 0) &
                call nc_put_att_str(this%db, trim(vars(j)%name), 'units', trim(vars(j)%units))
        end if
        ! Set statistic in CF-metadata convention        
        ! time cell_methods (mean or poitn) for this variable
        call nc_put_att_str(this%db, trim(vars(j)%name), 'cell_methods', trim(cm))

        if (vars(j)%has_min) then
            call nc_put_att_real(this%db, trim(vars(j)%name), 'valid_min', vars(j)%valid_min)
        end if
        if (vars(j)%has_max) then
            call nc_put_att_real(this%db, trim(vars(j)%name), 'valid_max', vars(j)%valid_max)
        end if
        if (vars(j)%has_missing) then
            call nc_put_att_real(this%db, trim(vars(j)%name), 'missing_value', vars(j)%missing_value)            
            !call nc_put_att_real(this%db, trim(vars(j)%name), '_FillValue', vars(j)%missing_value)
        end if
    end do    

    ! --- Attributes 
    if (present(title)) call nc_put_att_str(this%db, '', 'title', trim(title), global=.true.)
    call nc_put_att_str(this%db, '', 'Conventions', 'CF-1.10', global=.true.)

    call nc_put_att_str(this%db, 'time', 'units',    this%time_units)
    call nc_put_att_str(this%db, 'time', 'calendar', this%calendar_name)

    call nc_put_att_str(this%db, 'depth',   'standard_name',   'depth')
    call nc_put_att_str(this%db, 'depth',   'units',   'm')
    call nc_put_att_str(this%db, 'depth',   'positive','down')
    call nc_put_att_str(this%db, 'depth',   'axis',    'Z')
    call nc_put_att_str(this%db, 'depth_interface',   'long_name',   'depth of layer interfaces')
    call nc_put_att_str(this%db, 'depth_interface', 'units',   'm')
    call nc_put_att_str(this%db, 'depth_interface', 'positive','down')
    call nc_put_att_str(this%db, 'depth_interface', 'axis',    'Z')

    call nc_enddef(this%db)


    allocate(mask_out(this%N))
    mask_out(:) = 0.0_rk

    if (nsed_local > 0) then
        ! In file order: surface->bottom, sediments are last n_sed elements
        mask_out(this%N - nsed_local + 1 : this%N) = 1.0_rk
    end if

    call nc_write_real(this%db, 'sediment_mask', start=[1], count=[int(this%N,kind=4)], data_array=mask_out)
    deallocate(mask_out)


    ! --- Depth values (surface to bottom)
    allocate(z_out(this%N), zw_out(this%Ni))
    call flip_centers_b2s_to_s2b(grid%z,  z_out)         ! z(1..N) 
    call flip_interfaces_b2s_to_s2b(grid%z_w, zw_out)    ! z_w(0..N)

    call nc_write_real(this%db, 'depth',   start=[1], count=[int(this%N,kind=4)],   data_array=z_out)
    call nc_write_real(this%db, 'depth_interface', start=[1], count=[int(this%Ni,kind=4)],  data_array=zw_out)
    deallocate(z_out, zw_out)
  end subroutine outputwriter_open_file

  !========================
  ! Append one record per closed window
  !========================
  subroutine outputwriter_append_record(this, elapsed_time_s, window_len_s, is_mean, &
                                        centre_data, iface_data, scalar_data)
    class(OutputWriter), intent(inout) :: this
    integer(lk),         intent(in)    :: elapsed_time_s   ! seconds since sim start (window start)
    integer(lk),         intent(in)    :: window_len_s       ! seconds (window length)
    logical,             intent(in)    :: is_mean            ! .true. if statistic='mean'
    real(rk),            intent(in)    :: centre_data(:,:)   ! Data at layers' centres
    real(rk),            intent(in)    :: iface_data(:,:)    ! Data at interfaces
    real(rk),            intent(in)    :: scalar_data(:)     ! Scalars

    integer :: tidx, ivar
    integer, dimension(2) :: start2d, count2d
    integer, dimension(1) :: start1d, count1d
    real(rk) :: t_coord


    ! Timestamp
    t_coord = real(elapsed_time_s, rk)

    ! Increment time index
    this%time_index = this%time_index + 1_lk
    tidx = int(this%time_index, kind=4)

    call nc_write_real(this%db, 'time', start=[tidx],   count=[1], &
                       data_array=[t_coord])

    if (this%n_centre > 0) then
        if (size(centre_data, 1) /= this%N .or. size(centre_data, 2) /= this%n_centre) then
            error stop 'outputwriter_append_record: centre_data has wrong shape.'
        end if
    end if

    if (this%n_iface > 0) then
        if (size(iface_data, 1) /= this%Ni .or. size(iface_data, 2) /= this%n_iface) then
            error stop 'outputwriter_append_record: iface_data has wrong shape.'
        end if
    end if

    if (this%n_scalar > 0) then
        if (size(scalar_data) /= this%n_scalar) then
            error stop 'outputwriter_append_record: scalar_data has wrong length.'
        end if
    end if


    ! Data at layers' centres
    if (this%n_centre > 0) then
      start2d = [tidx, 1]
      count2d = [1, int(this%N, kind=4)]
      
      do ivar = 1, this%n_centre  
        ! Flip from bottom-surface to surface-bottom
        call flip_centers_b2s_to_s2b(centre_data(:, ivar), this%center_row)      

        call nc_write_real(this%db, this%centre_varids(ivar)%name, start=start2d, count=count2d, &
                          data_array=reshape(this%center_row, [1, int(this%N,kind=4)]))
      end do
    end if 
    ! For data at the interfaces
    if (this%n_iface > 0) then
      start2d = [tidx, 1]
      count2d = [1, int(this%Ni, kind=4)]

      do ivar = 1, this%n_iface
          call flip_interfaces_b2s_to_s2b(iface_data(:, ivar), this%iface_row)       

          call nc_write_real(this%db, this%iface_varids(ivar)%name, start=start2d, count=count2d, &
                            data_array=reshape(this%iface_row, [1, int(this%Ni,kind=4)]))
      end do
    end if

    if (this%n_scalar > 0) then
      start1d = [tidx]
      count1d = [1]

      do ivar = 1, this%n_scalar
          call nc_write_real(this%db, this%scalar_varids(ivar)%name, &
                            start=start1d, count=count1d,                 &
                            data_array=[scalar_data(ivar)])
      end do 
    end if   
  end subroutine outputwriter_append_record

  !========================
  ! Sync (flush) to disk
  !========================
  subroutine outputwriter_sync_file(this)
    class(OutputWriter), intent(inout) :: this
    call nc_sync(this%db)
  end subroutine outputwriter_sync_file

  !========================
  ! Close
  !========================
  subroutine outputwriter_close_file(this, do_sync)
    class(OutputWriter), intent(inout) :: this
    logical, intent(in), optional      :: do_sync
    if (present(do_sync)) then
      if (do_sync) call nc_sync(this%db)
    end if
    call nc_close(this%db)
  end subroutine outputwriter_close_file


  pure function cf_cell_methods_from_stat(statistic) result(cm)
    character(*), intent(in) :: statistic       ! 'mean' or 'instant'
    character(len=:), allocatable :: cm
    select case (trim(adjustl(statistic)))
    case ('mean');     cm = 'time: mean'
    case ('instant');  cm = 'time: point'
    case default
       cm = 'time: point'   ! safe default
    end select
  end function cf_cell_methods_from_stat

  !========================
  ! Flip helpers (physics: bottom→surface; file: surface→bottom)
  !========================

  pure subroutine flip_centers_b2s_to_s2b(src, dst)
    real(rk), intent(in)  :: src(:)    ! size N (bottom→surface)
    real(rk), intent(out) :: dst(:)    ! size N (surface→bottom)
    integer :: N, k
    N = size(src); if (size(dst) /= N) error stop 'flip_centers: size mismatch'
    do k = 1, N
      dst(k) = src(N - k + 1)
    end do
  end subroutine flip_centers_b2s_to_s2b

  pure subroutine flip_interfaces_b2s_to_s2b(src, dst)
    ! Accepts src indexed either 0..N or 1..N+1, both bottom to surface
    real(rk), intent(in)  :: src(:)    ! size N+1 (bottom to surface)
    real(rk), intent(out) :: dst(:)    ! size N+1 (surface to bottom)
    integer :: Np1, lb, k
    Np1 = size(src); if (size(dst) /= Np1) error stop 'flip_interfaces: size mismatch'
    lb  = lbound(src,1)
    select case (lb)
    case (0)
      ! src(0) bottom … src(N) surface  to dst( N+1: bottom … 1:surface )
      do k = 1, Np1
        dst(k) = src(Np1 - k)   ! since last index is N = Np1-1
      end do
    case (1)
      ! src(1) bottom … src(N+1) surface
      do k = 1, Np1
        dst(k) = src(Np1 - k + 1)
      end do
    case default
      error stop 'flip_interfaces: unexpected lower bound'
    end select
  end subroutine flip_interfaces_b2s_to_s2b

end module output_writer
