! src/io/output_writer.F90
module output_writer
  use iso_fortran_env, only: error_unit
  use precision_types,  only: rk, lk
  use netcdf,           only: NF90_UNLIMITED
  use netcdf_io,        only: NcFile, nc_create, nc_redef, nc_enddef, nc_close, nc_sync, &
                              nc_def_dim, nc_def_var_real, nc_def_var_double, &
                              nc_put_att_str, nc_put_att_real, nc_write_real
  use grids,            only: VerticalGrid
  implicit none
  private
  public :: OutputWriter, outputwriter_open_file, outputwriter_append_record, outputwriter_close_file
  public :: cf_cell_methods_from_stat


  !========================
  ! Writer type
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

     ! Biogeochem tracers
     integer                   :: n_bio = 0
     character(:), allocatable :: bio_names(:)
     real(rk), allocatable :: bio_row(:,:)      ! (N, n_bio)

     ! working row buffers in file order (surface→bottom), allocated once
     real(rk), allocatable :: temp_row(:)   ! size N
     real(rk), allocatable :: kz_row(:)     ! size N+1
  contains
     procedure :: open_file    => outputwriter_open_file
     procedure :: append_record=> outputwriter_append_record
     procedure :: close_file   => outputwriter_close_file
  end type OutputWriter

contains  

  !========================
  ! Open & define the file (write coords & attrs)
  !========================
  subroutine outputwriter_open_file(this, path, grid_phys, interval_statistic, time_units, calendar_name, title,n_bio, bio_names)
    class(OutputWriter), intent(inout) :: this
    character(*),        intent(in)    :: path
    type(VerticalGrid),  intent(in)    :: grid_phys        ! bottom→surface order internally
    character(*),        intent(in)    :: interval_statistic
    character(*),        intent(in)    :: time_units       ! "seconds since ..."
    character(*),        intent(in)    :: calendar_name
    character(*),        intent(in), optional :: title
    integer,             intent(in), optional :: n_bio
    character(*),        intent(in), optional :: bio_names(:)

    integer(lk) :: nb
    integer :: dim_t, dim_z, dim_zw
    integer :: var_time, var_z, var_zw, var_temp, var_kz, var_bio 
    integer :: j
    character(len=:), allocatable :: cm
    real(rk), allocatable :: z_out(:), zw_out(:)

    ! sizes
    this%N  = grid_phys%nz
    this%Ni = grid_phys%nz + 1

    this%n_bio = 0_lk
    if (present(n_bio)) this%n_bio = n_bio

    if (present(bio_names)) then
        if (this%n_bio > 0_lk .and. this%n_bio /= size(bio_names, kind=lk)) then
            error stop 'outputwriter_open_file: n_bio does not match size(bio_names).'
        end if
        allocate(this%bio_names, source=bio_names)   ! ✔ length + shape taken from dummy
        this%n_bio = size(bio_names, kind=lk)        ! in case n_bio was not passed
    end if

    allocate(this%temp_row(this%N), this%kz_row(this%Ni))
    if (this%n_bio > 0_lk) allocate(this%bio_row(this%N, this%n_bio))

    this%filename           = trim(path)
    this%interval_statistic = trim(interval_statistic)
    this%time_units         = trim(time_units)
    this%calendar_name      = trim(calendar_name)
    this%time_index         = 0_lk


    ! Create & define structure
    call nc_create(this%db, this%filename, overwrite=.true.)
    call nc_redef(this%db)

    dim_t  = nc_def_dim(this%db, 'time',  NF90_UNLIMITED)
    dim_z  = nc_def_dim(this%db, 'z',     int(this%N, kind=4))
    dim_zw = nc_def_dim(this%db, 'z_w',   int(this%Ni, kind=4))

    var_time = nc_def_var_double(this%db, 'time', [dim_t])              ! NOTE: uses rk; prefer rk=real64
    var_z    = nc_def_var_real(this%db, 'z',    [dim_z])
    var_zw   = nc_def_var_real(this%db, 'z_w',  [dim_zw])

    var_temp = nc_def_var_real(this%db, 'temp', [dim_t, dim_z])
    var_kz   = nc_def_var_real(this%db, 'Kz',   [dim_t, dim_zw])

    ! Biogeochemical state variables: one var per tracer
    if (this%n_bio > 0_lk) then
       do j = 1, int(this%n_bio, kind=4)
          var_bio = nc_def_var_real(this%db, trim(this%bio_names(j)), [dim_t, dim_z])
          ! For now, no units/long_name; you’ll add them later.
       end do
    end if

    ! --- Attributes (still in define mode)
    if (present(title)) call nc_put_att_str(this%db, '', 'title', trim(title), global=.true.)
    call nc_put_att_str(this%db, '', 'Conventions', 'CF-1.10', global=.true.)

    call nc_put_att_str(this%db, 'time', 'units',    this%time_units)
    call nc_put_att_str(this%db, 'time', 'calendar', this%calendar_name)

    call nc_put_att_str(this%db, 'z',   'standard_name',   'depth')
    call nc_put_att_str(this%db, 'z',   'units',   'm')
    call nc_put_att_str(this%db, 'z',   'positive','down')
    call nc_put_att_str(this%db, 'z',   'axis',    'Z')
    call nc_put_att_str(this%db, 'z_w',   'standard_name',   'depth of layer interfaces')
    call nc_put_att_str(this%db, 'z_w', 'units',   'm')
    call nc_put_att_str(this%db, 'z_w', 'positive','down')
    call nc_put_att_str(this%db, 'z_w', 'axis',    'Z')

    call nc_put_att_str(this%db, 'temp', 'long_name', 'potential temperature')
    call nc_put_att_str(this%db, 'temp', 'units',     'degree_Celsius')
    cm = cf_cell_methods_from_stat(this%interval_statistic)
    call nc_put_att_str(this%db, 'temp', 'cell_methods', cm)

    call nc_put_att_str(this%db, 'Kz', 'long_name', 'vertical diffusivity')
    call nc_put_att_str(this%db, 'Kz', 'units',     'm2 s-1')
    call nc_put_att_str(this%db, 'Kz', 'cell_methods', cm)

    call nc_enddef(this%db)

    ! --- Coordinate values (surface→bottom)
    allocate(z_out(this%N), zw_out(this%Ni))
    call flip_centers_b2s_to_s2b(grid_phys%z,  z_out)         ! z(1..N) bottom→surface → z_out(1..N) surface→bottom
    call flip_interfaces_b2s_to_s2b(grid_phys%z_w, zw_out)    ! z_w(0..N or 1..N+1) → zw_out(1..N+1) surface→bottom

    call nc_write_real(this%db, 'z',   start=[1], count=[int(this%N,kind=4)],   data_array=z_out)
    call nc_write_real(this%db, 'z_w', start=[1], count=[int(this%Ni,kind=4)],  data_array=zw_out)
  end subroutine outputwriter_open_file

  !========================
  ! Append one record per closed window
  !========================
  subroutine outputwriter_append_record(this, t_window_start_s, window_len_s, temp_phys, kz_phys, is_mean, bio_phys)
    class(OutputWriter), intent(inout) :: this
    integer(lk),         intent(in)    :: t_window_start_s   ! seconds since sim start (window start)
    integer(lk),         intent(in)    :: window_len_s       ! seconds (window length)
    real(rk),            intent(in)    :: temp_phys(:)       ! N, bottom→surface
    real(rk),            intent(in)    :: kz_phys(:)         ! N+1, bottom→surface (lb 0 or 1)
    logical,             intent(in)    :: is_mean            ! .true. if statistic='mean'
    real(rk),            intent(in), optional :: bio_phys(:,:)   ! (N, n_bio), bottom to surface

    integer :: tidx, j
    real(rk) :: t_coord

    if (size(temp_phys) /= this%N)  error stop 'append_record: temp size mismatch'
    if (size(kz_phys)   /= this%Ni) error stop 'append_record: Kz size mismatch'

    if (present(bio_phys) .and. this%n_bio > 0_lk) then
       if (size(bio_phys,1) /= this%N .or. size(bio_phys,2) /= this%n_bio) then
          error stop 'append_record: bio_phys size mismatch'
       end if
    end if

    ! 1) Flip rows into file order (surface→bottom)
    call flip_centers_b2s_to_s2b(temp_phys, this%temp_row)
    call flip_interfaces_b2s_to_s2b(kz_phys, this%kz_row)

    if (present(bio_phys) .and. this%n_bio > 0_lk) then
       do j = 1, int(this%n_bio, kind=4)
          call flip_centers_b2s_to_s2b(bio_phys(:, j), this%bio_row(:, j))
       end do
    end if

    ! 2) Time stamp...
    if (is_mean) then
       t_coord = real(t_window_start_s + window_len_s/2_lk, rk)
    else
       t_coord = real(t_window_start_s + window_len_s, rk)
    end if

    ! 3) Write one slice at next time index
    this%time_index = this%time_index + 1_lk
    tidx = int(this%time_index, kind=4)

    call nc_write_real(this%db, 'time', start=[tidx],   count=[1], &
                       data_array=[t_coord])
    call nc_write_real(this%db, 'temp', start=[tidx,1], count=[1, int(this%N,kind=4)], &
                       data_array=reshape(this%temp_row, [1, int(this%N,kind=4)]))
    call nc_write_real(this%db, 'Kz',   start=[tidx,1], count=[1, int(this%Ni,kind=4)], &
                       data_array=reshape(this%kz_row,   [1, int(this%Ni,kind=4)]))

    if (this%n_bio > 0_lk .and. present(bio_phys)) then
       do j = 1, int(this%n_bio, kind=4)
          call nc_write_real(this%db, trim(this%bio_names(j)), &
                             start=[tidx,1], count=[1, int(this%N,kind=4)], &
                             data_array=reshape(this%bio_row(:,j), [1, int(this%N,kind=4)]))
       end do
    end if
  end subroutine outputwriter_append_record

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
    ! Accepts src indexed either 0..N or 1..N+1, both bottom→surface
    real(rk), intent(in)  :: src(:)    ! size N+1 (bottom→surface)
    real(rk), intent(out) :: dst(:)    ! size N+1 (surface→bottom)
    integer :: Np1, lb, k
    Np1 = size(src); if (size(dst) /= Np1) error stop 'flip_interfaces: size mismatch'
    lb  = lbound(src,1)
    select case (lb)
    case (0)
      ! src(0) bottom … src(N) surface  → dst( N+1→bottom … 1→surface )
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
