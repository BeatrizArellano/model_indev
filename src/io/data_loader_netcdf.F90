module data_loader_netcdf
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite   
   
   use cf_time_utils,   only: parse_cf_time, seconds_since_datetime_file
   use data_types,      only: DataSpec, DataVarSeries, DATA_INPUT_FILE
   use find_utils,      only: find_name, has_name, argmin_abs_vec
   use geo_utils,       only: LocationInfo, simple_distance_deg
   use netcdf
   use netcdf_io,       only: NcFile, nc_open, nc_close, nc_has_var, nc_check, &
                              nc_var_dims, nc_read_real_1d, nc_get_att_str   
   use precision_types, only: rk, lk
   use str_utils,       only: inttostr, realtostr, list_to_str, append_string   
   use time_types,      only: CFUnits, TimeAxis, CFCalendar, DateTime
   use time_utils,      only: detect_frequency, index_at_or_before, index_at_or_after, check_time_monotonic

   implicit none
   private

   public :: NetcdfScan
   public :: scan_netcdf_file
   public :: build_netcdf_year_windows
   public :: load_netcdf_series
   public :: read_netcdf_timeseries_at_point

   integer(lk), parameter :: INF_EDGE = huge(1_lk)

   type :: NetcdfScan
      character(:), allocatable :: path

      character(:), allocatable :: time_name
      character(:), allocatable :: time_dim

      character(:), allocatable :: lat_name
      character(:), allocatable :: lat_dim
      character(:), allocatable :: lon_name
      character(:), allocatable :: lon_dim

      logical :: has_latlon = .false.
      logical :: is_point   = .true.
      integer :: yi = 1
      integer :: xi = 1

      real(rk), allocatable :: lat(:)
      real(rk), allocatable :: lon(:)
      character(:), allocatable :: lon_convention
      real(rk) :: lat_found = 0.0_rk
      real(rk) :: lon_found = 0.0_rk

      type(CFUnits)    :: u
      type(CFCalendar) :: cal
      type(TimeAxis)   :: axis

      integer :: i0 = 1
      integer :: i1 = 1
      real(rk) :: median_dt = 0.0_rk
      logical  :: is_regular = .true.
      real(rk) :: rel_max_dev = 0.0_rk
      integer(lk) :: sim_offset = 0_lk

      character(:), allocatable :: vars_present(:)
      character(:), allocatable :: vars_missing(:)

      integer :: y0 = 0, mon0 = 0, d0 = 0, h0 = 0, mi0 = 0, s0 = 0
      integer :: y1 = 0, mon1 = 0, d1 = 0, h1 = 0, mi1 = 0, s1 = 0
   end type NetcdfScan

contains

   subroutine scan_netcdf_file(path, specs, location, start_datetime, end_datetime, &
                               calendar_default, scan, ok, errmsg, max_sep_deg)
      character(*),        intent(in)    :: path
      type(DataSpec),      intent(inout) :: specs(:)
      type(LocationInfo),  intent(in)    :: location
      type(DateTime),      intent(in)    :: start_datetime, end_datetime
      character(*),        intent(in)    :: calendar_default
      type(NetcdfScan),    intent(out)   :: scan
      logical,             intent(out)   :: ok
      character(*),        intent(out)   :: errmsg
      real(rk), optional,  intent(in)    :: max_sep_deg

      type(NcFile) :: db
      character(:), allocatable :: required_vars(:)
      character(:), allocatable :: time_name
      character(:), allocatable :: lat_name, lon_name
      character(:), allocatable :: units_attr, cal_attr
      character(:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)
      real(rk), allocatable :: time_orig(:)
      logical :: pres, ok_parse, lok
      logical :: has_time, has_lat, has_lon
      logical :: is_monotonic, has_equal_consecutive
      integer :: ntime, i_notmon, i
      real(rk) :: t_start, t_end, tol_cover

      ok = .false.
      errmsg = ''
      call clear_scan(scan)
      scan%path = trim(path)

      call collect_required_vars(specs, path, required_vars)

      call nc_open(db, path)

      call find_time_name(db, specs, path, time_name, has_time)
      if (.not. has_time) then
         errmsg = 'NetCDF file '//trim(path)//' is missing the required time coordinate.'
         call nc_close(db)
         return
      end if
      scan%time_name = trim(time_name)

      call nc_get_att_str(db, scan%time_name, 'units', units_attr, pres)
      if (.not. pres) then
         errmsg = 'CF time: missing units attribute on variable '//trim(scan%time_name)//' in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      call nc_get_att_str(db, scan%time_name, 'calendar', cal_attr, pres)
      if (.not. pres) cal_attr = ''

      call parse_cf_time(units_attr, cal_attr, trim(calendar_default), scan%u, scan%cal, ok_parse, errmsg)
      if (.not. ok_parse) then
         call nc_close(db)
         return
      end if

      call nc_var_dims(db, scan%time_name, dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'CF time variable '//trim(scan%time_name)//' must be 1-D in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      scan%time_dim = trim(dimnames(1))

      ntime = dimlens(1)
      if (ntime < 1) then
         errmsg = 'CF time variable '//trim(scan%time_name)//' is empty in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      allocate(time_orig(ntime))
      call nc_read_real_1d(db, scan%time_name, time_orig)

      scan%axis%cal = scan%cal
      scan%axis%u   = scan%u
      allocate(scan%axis%t_s(ntime))
      scan%axis%t_s     = time_orig * scan%u%timeunit_to_seconds
      scan%axis%t_first = scan%axis%t_s(1)
      scan%axis%t_last  = scan%axis%t_s(ntime)

      if (any(.not. ieee_is_finite(scan%axis%t_s))) then
         errmsg = 'The time coordinate contains NaN/Inf values in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      call check_time_monotonic(scan%axis%t_s, is_monotonic, has_equal_consecutive, i_notmon)
      if (.not. is_monotonic) then
         errmsg = 'CF time is not sorted in '//trim(path)//': first inversion at index '// &
                  trim(adjustl(inttostr(i_notmon)))//'-'//trim(adjustl(inttostr(i_notmon + 1)))//'.'
         call nc_close(db)
         return
      end if

      if (has_equal_consecutive) then
         errmsg = 'CF time contains duplicate timestamps in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      call detect_frequency(scan%axis%t_s, scan%median_dt, scan%is_regular, scan%rel_max_dev)
      if (scan%median_dt <= 0.0_rk) then
         errmsg = 'Detected non-positive NetCDF time step in '//trim(path)//'.'
         call nc_close(db)
         return
      end if

      call find_horizontal_coordinates(db, location, scan, lok, errmsg, max_sep_deg)
      if (.not. lok) then
         ok = .false.
         call nc_close(db)
         return
      end if

      scan%y0   = start_datetime%year
      scan%mon0 = start_datetime%month
      scan%d0   = start_datetime%day
      scan%h0   = start_datetime%hour
      scan%mi0  = start_datetime%minute
      scan%s0   = start_datetime%second

      scan%y1   = end_datetime%year
      scan%mon1 = end_datetime%month
      scan%d1   = end_datetime%day
      scan%h1   = end_datetime%hour
      scan%mi1  = end_datetime%minute
      scan%s1   = end_datetime%second

      t_start = seconds_since_datetime_file(scan%cal, scan%u, scan%y0, scan%mon0, scan%d0, scan%h0, scan%mi0, scan%s0)
      t_end   = seconds_since_datetime_file(scan%cal, scan%u, scan%y1, scan%mon1, scan%d1, scan%h1, scan%mi1, scan%s1)

      tol_cover = 0.999_rk * scan%median_dt
      if (t_start < scan%axis%t_first - tol_cover .or. t_end > scan%axis%t_last + tol_cover .or. t_end < t_start) then
         errmsg = 'NetCDF file '//trim(path)//' does not cover the requested simulation period.'
         call nc_close(db)
         return
      end if

      scan%sim_offset = nint(t_start, kind=lk)
      scan%i0 = max(1, index_at_or_before(scan%axis%t_s, t_start))
      scan%i1 = max(scan%i0, index_at_or_before(scan%axis%t_s, t_end))

      call clear_char_list(scan%vars_present)
      call clear_char_list(scan%vars_missing)

      ! Preserve the existing scalar pathway exactly as before. Profile variables
      ! are validated separately because their vertical coordinate is spec-specific.
      do i = 1, size(required_vars)
         call check_var_dims_cf(db, trim(required_vars(i)), scan%has_latlon, scan%time_dim, &
                                scan%lat_dim, scan%lon_dim, scan%vars_present, scan%vars_missing)
      end do

      if (allocated(scan%vars_missing) .and. size(scan%vars_missing) > 0) then
         errmsg = 'Missing or incompatible NetCDF variables in '//trim(path)//': '// &
                  trim(list_to_str(scan%vars_missing, ', '))
         call nc_close(db)
         return
      end if

      do i = 1, size(required_vars)
         if (.not. is_var_valid_in_period(db, trim(required_vars(i)), scan%time_dim, &
                                          scan%lat_dim, scan%lon_dim, scan%has_latlon, &
                                          scan%i0, scan%i1, scan%yi, scan%xi)) then
            errmsg = 'Variable '//trim(required_vars(i))// &
                     ' contains NaN/Inf values in '//trim(path)// &
                     ' during the simulation period. Verify the forcing data for this location. '
            
            call nc_close(db)
            ok = .false. 
            return
         end if
      end do

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
         if (.not. specs(i)%is_profile) cycle
         if (.not. allocated(specs(i)%path)) cycle
         if (trim(specs(i)%path) /= trim(path)) cycle

         call prepare_profile_spec(db, specs(i), scan%time_dim, scan%has_latlon, &
                                   scan%lat_dim, scan%lon_dim, lok, errmsg)
         if (.not. lok) then
            call nc_close(db)
            ok = .false.
            return
         end if

         call append_string(scan%vars_present, trim(specs(i)%source_var))
      end do

      ok = .true.
      call nc_close(db)
   end subroutine scan_netcdf_file


   subroutine build_netcdf_year_windows(spec, scan, source_y_start, source_y_end, &
                                        window_start_datetime, window_end_datetime)
      type(DataSpec),   intent(inout) :: spec
      type(NetcdfScan), intent(in)    :: scan
      integer,          intent(in)    :: source_y_start, source_y_end
      type(DateTime),   intent(in)    :: window_start_datetime, window_end_datetime

      integer  :: ny, k, y
      real(rk) :: t0, t1, t1_exclusive, tol_cover
      real(rk) :: tw_start, tw_end

      ny = source_y_end - source_y_start + 1
      if (ny < 1) error stop 'build_netcdf_year_windows: invalid source year range.'

      if (allocated(spec%idx_window)) deallocate(spec%idx_window)
      allocate(spec%idx_window(2, ny))

      tol_cover = 0.999_rk * scan%median_dt

      tw_start = seconds_since_datetime_file(scan%cal, scan%u, &
                                       window_start_datetime%year, window_start_datetime%month, window_start_datetime%day, &
                                       window_start_datetime%hour, window_start_datetime%minute, window_start_datetime%second)

      tw_end = seconds_since_datetime_file(scan%cal, scan%u, &
                                          window_end_datetime%year, window_end_datetime%month, window_end_datetime%day, &
                                          window_end_datetime%hour, window_end_datetime%minute, window_end_datetime%second)

      do k = 1, ny
         y = source_y_start + k - 1

         t0 = seconds_since_datetime_file(scan%cal, scan%u, y,     1, 1, 0, 0, 0)
         t1 = seconds_since_datetime_file(scan%cal, scan%u, y + 1, 1, 1, 0, 0, 0)

         t0 = max(t0, tw_start)
         t1 = min(t1, tw_end)
         
         if (t0 < scan%axis%t_first - tol_cover .or. &
            t1 > scan%axis%t_last  + tol_cover) then
            error stop 'build_netcdf_year_windows: requested interval is not covered by NetCDF file.'
         end if

         ! Use a half-open interval [t0, t1), so that Jan 1 of the next year
         ! is not included in the current year window if it exists in the file.
         t1_exclusive = t1 - 1.0_rk

         spec%idx_window(1, k) = max(1, index_at_or_after(scan%axis%t_s, t0))
         spec%idx_window(2, k) = max(spec%idx_window(1, k), index_at_or_before(scan%axis%t_s, t1_exclusive))
      end do
   end subroutine build_netcdf_year_windows


   subroutine load_netcdf_series(series, scan, db, spec, i0, i1)
      type(DataVarSeries), intent(inout) :: series
      type(NetcdfScan),    intent(in)    :: scan
      type(NcFile),        intent(in)    :: db
      type(DataSpec),      intent(in)    :: spec
      integer,             intent(in)    :: i0, i1

      integer :: nt, i
      integer(lk) :: dt_last

      nt = max(0, i1 - i0 + 1)
      if (nt <= 0) error stop 'load_netcdf_series: empty time window for '//trim(spec%name)

      if (allocated(series%t_axis)) deallocate(series%t_axis)
      if (allocated(series%t_edge)) deallocate(series%t_edge)
      if (allocated(series%values)) deallocate(series%values)

      allocate(series%t_axis(nt))
      allocate(series%values(nt))

      series%name     = trim(spec%name)
      series%units    = trim(spec%units)
      series%is_const = .false.
      series%n        = nt
      series%idx      = 1

      series%t_axis = nint(scan%axis%t_s(i0:i1), kind=lk)

      call read_netcdf_timeseries_at_point(db, trim(spec%source_var), trim(scan%time_dim), i0, i1, &
                                           scan%has_latlon, trim(scan%lat_dim), trim(scan%lon_dim), &
                                           scan%yi, scan%xi, series%values)

      allocate(series%t_edge(nt + 1))

      ! Left or step-ahead convention:
      ! value(i) is valid over [t_axis(i), t_axis(i+1)).
      ! The timestamp marks the beginning of the interval, not its midpoint.
      if (nt >= 2) then
         do i = 1, nt
            series%t_edge(i) = series%t_axis(i)
         end do

         dt_last = max(1_lk, series%t_axis(nt) - series%t_axis(nt - 1))
         series%t_edge(nt + 1) = series%t_axis(nt) + dt_last

         series%t_next = series%t_edge(2)
      else
         dt_last = max(1_lk, nint(scan%median_dt, kind=lk))

         series%t_edge(1) = series%t_axis(1)
         series%t_edge(2) = series%t_axis(1) + dt_last

         series%t_next = INF_EDGE
      end if
   end subroutine load_netcdf_series


   subroutine read_netcdf_timeseries_at_point(db, varname, time_dim, i0, i1, &
                                              has_latlon, lat_dim, lon_dim, yi, xi, out)
      type(NcFile), intent(in)  :: db
      character(*), intent(in)  :: varname, time_dim, lat_dim, lon_dim
      logical,      intent(in)  :: has_latlon
      integer,      intent(in)  :: i0, i1, yi, xi
      real(rk),     intent(out) :: out(:)

      integer :: vid, ndims, dimids(NF90_MAX_VAR_DIMS), xtype, natts
      character(:), allocatable :: dnames(:)
      integer, allocatable :: dlens(:)
      integer :: itime, ilat, ilon, nt
      integer :: start(NF90_MAX_VAR_DIMS), count(NF90_MAX_VAR_DIMS)
      logical :: uses_latlon, ok_dims
      character(len=512) :: errmsg_dims

      call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
      call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                    'inquire_variable('//trim(varname)//')')

      if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'read_netcdf_timeseries_at_point: '//trim(varname)//' is character')

      call nc_var_dims(db, varname, dnames, dlens)

      call resolve_timeseries_dims(db, varname, time_dim, has_latlon, lat_dim, lon_dim, &
                                   uses_latlon, itime, ilat, ilon, ok_dims, errmsg_dims)
      if (.not. ok_dims) then
         write(*,'(A)') trim(errmsg_dims)
         call nc_check(NF90_EBADDIM, &
                        'read_netcdf_timeseries_at_point: '//trim(varname))
      end if

      nt = i1 - i0 + 1
      if (nt < 1) call nc_check(NF90_EEDGE, 'empty time window for '//trim(varname))
      if (i0 < 1 .or. i1 > dlens(itime)) call nc_check(NF90_EEDGE, 'time indices out of range for '//trim(varname))
      if (size(out) /= nt) call nc_check(NF90_EEDGE, 'output size mismatch for '//trim(varname))

      start(1:ndims) = 1
      count(1:ndims) = 1

      start(itime) = i0
      count(itime) = nt

      if (uses_latlon) then
         if (yi < 1 .or. yi > dlens(ilat) .or. &
            xi < 1 .or. xi > dlens(ilon)) then
            call nc_check(NF90_EEDGE, &
                           'selected lat/lon index out of range for '//trim(varname))
         end if

         start(ilat) = yi
         start(ilon) = xi
      end if

      call nc_check(nf90_get_var(db%ncid, vid, out, start=start(1:ndims), count=count(1:ndims)), &
                    'get_var slice '//trim(varname))
   end subroutine read_netcdf_timeseries_at_point


   subroutine collect_required_vars(specs, path, names)
      type(DataSpec), intent(in) :: specs(:)
      character(*),   intent(in) :: path
      character(:), allocatable, intent(out) :: names(:)

      integer :: i, n, maxlen
      logical :: found

      maxlen = 1
      n = 0

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
         if (specs(i)%is_profile) cycle
         if (.not. allocated(specs(i)%path)) cycle
         if (trim(specs(i)%path) /= trim(path)) cycle
         if (.not. allocated(specs(i)%source_var)) cycle
         maxlen = max(maxlen, len_trim(specs(i)%source_var))
      end do

      allocate(character(len=maxlen) :: names(size(specs)))
      names = ''

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
         if (specs(i)%is_profile) cycle
         if (.not. allocated(specs(i)%path)) cycle
         if (trim(specs(i)%path) /= trim(path)) cycle

         found = .false.
         if (n > 0) found = has_name(names(1:n), trim(specs(i)%source_var))

         if (.not. found) then
            n = n + 1
            names(n) = trim(specs(i)%source_var)
         end if
      end do

      if (n == 0) then
         deallocate(names)
         allocate(character(len=1) :: names(0))
      else if (n < size(names)) then
         names = names(1:n)
      end if
   end subroutine collect_required_vars


   subroutine find_time_name(db, specs, path, time_name, found)
      type(NcFile),   intent(in)  :: db
      type(DataSpec), intent(in)  :: specs(:)
      character(*),   intent(in)  :: path
      character(:), allocatable, intent(out) :: time_name
      logical,        intent(out) :: found

      integer :: i

      time_name = ''
      found = .false.

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
         if (.not. allocated(specs(i)%path)) cycle
         if (trim(specs(i)%path) /= trim(path)) cycle
         if (allocated(specs(i)%time_var) .and. len_trim(specs(i)%time_var) > 0) then
            time_name = trim(specs(i)%time_var)
            found = nc_has_var(db, time_name)
            return
         end if
      end do

      if (nc_has_var(db, 'time')) then
         time_name = 'time'
         found = .true.
      end if
   end subroutine find_time_name


   subroutine find_horizontal_coordinates(db, location, scan, ok, errmsg, max_sep_deg)
      type(NcFile),       intent(in)    :: db
      type(LocationInfo), intent(in)    :: location
      type(NetcdfScan),   intent(inout) :: scan
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg
      real(rk), optional, intent(in)    :: max_sep_deg

      character(:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)
      logical :: has_lat, has_lon
      real(rk) :: lonq, sep_limit, sep_deg
      real(rk), parameter :: default_max_sep_deg = 0.25_rk

      ok = .false.
      errmsg = ''

      has_lat = nc_has_var(db, 'latitude') .or. nc_has_var(db, 'lat')
      has_lon = nc_has_var(db, 'longitude') .or. nc_has_var(db, 'lon') .or. nc_has_var(db, 'long')

      if (nc_has_var(db, 'latitude')) then
         scan%lat_name = 'latitude'
      else
         scan%lat_name = 'lat'
      end if

      if (nc_has_var(db, 'longitude')) then
         scan%lon_name = 'longitude'
      else if (nc_has_var(db, 'lon')) then
         scan%lon_name = 'lon'
      else
         scan%lon_name = 'long'
      end if

      scan%has_latlon = has_lat .and. has_lon
      scan%is_point = .not. scan%has_latlon

      if (.not. scan%has_latlon) then
         scan%lon_convention = '-'
         scan%lat_found = location%lat
         scan%lon_found = location%lon
         scan%lat_dim = ''
         scan%lon_dim = ''
         scan%yi = 1
         scan%xi = 1
         ok = .true.
         return
      end if

      call nc_var_dims(db, scan%lat_name, dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'Latitude coordinate must be 1-D; curvilinear grids are not supported yet.'
         return
      end if
      scan%lat_dim = trim(dimnames(1))
      allocate(scan%lat(dimlens(1)))
      call nc_read_real_1d(db, scan%lat_name, scan%lat)

      call nc_var_dims(db, scan%lon_name, dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'Longitude coordinate must be 1-D; curvilinear grids are not supported yet.'
         return
      end if
      scan%lon_dim = trim(dimnames(1))
      allocate(scan%lon(dimlens(1)))

      if (trim(scan%lat_dim) == trim(scan%lon_dim)) then
         errmsg = 'Latitude and longitude coordinates must use distinct dimensions; '// &
                  'paired horizontal coordinates are not supported yet.'
         return
      end if

      call nc_read_real_1d(db, scan%lon_name, scan%lon)

      if (any(.not. ieee_is_finite(scan%lat)) .or. any(.not. ieee_is_finite(scan%lon))) then
         errmsg = 'Latitude/longitude coordinates contain NaN/Inf values.'
         return
      end if

      if (minval(scan%lon) >= 0.0_rk .and. maxval(scan%lon) <= 360.0_rk) then
         scan%lon_convention = '0_360'
         lonq = modulo(location%lon, 360.0_rk)
         if (lonq < 0.0_rk) lonq = lonq + 360.0_rk
      else
         scan%lon_convention = '-180_180'
         lonq = modulo(location%lon + 180.0_rk, 360.0_rk) - 180.0_rk
      end if

      scan%yi = argmin_abs_vec(scan%lat - location%lat)
      scan%xi = argmin_abs_vec(scan%lon - lonq)
      scan%lat_found = scan%lat(scan%yi)
      scan%lon_found = scan%lon(scan%xi)

      if (present(max_sep_deg)) then
         sep_limit = max_sep_deg
      else
         sep_limit = default_max_sep_deg
      end if

      sep_deg = simple_distance_deg(location%lat, lonq, scan%lat_found, scan%lon_found)
      if (sep_deg > sep_limit) then
         errmsg = 'No nearby NetCDF grid point: nearest is ~'//realtostr(sep_deg, 3)// &
                  ' deg from requested site (> '//realtostr(sep_limit, 3)//' deg).'
         return
      end if

      ok = .true.
   end subroutine find_horizontal_coordinates


   subroutine resolve_timeseries_dims(db, vname, time_dim, has_latlon, lat_dim, lon_dim, &
                                   uses_latlon, itime, ilat, ilon, ok, errmsg)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname
      character(*), intent(in) :: time_dim
      logical,      intent(in) :: has_latlon
      character(*), intent(in) :: lat_dim, lon_dim

      logical,      intent(out) :: uses_latlon
      integer,      intent(out) :: itime, ilat, ilon
      logical,      intent(out) :: ok
      character(*), intent(out) :: errmsg

      character(:), allocatable :: dnames(:)
      integer, allocatable :: dlens(:)
      integer :: ndims

      ok = .false.
      errmsg = ''

      uses_latlon = .false.
      itime = 0
      ilat  = 0
      ilon  = 0

      call nc_var_dims(db, vname, dnames, dlens)

      ndims = size(dnames)

      itime = find_name(dnames, trim(time_dim))
      if (itime <= 0) then
         errmsg = 'Variable '//trim(vname)//' does not contain time dimension '// &
                  trim(time_dim)//'.'
         return
      end if

      ! Pure time-dependent variable: var(time)
      if (ndims == 1) then
         ok = .true.
         return
      end if

      ! Time-dependent spatial variable: var(time, lat, lon)
      ! Dimension order does not matter.
      if (ndims == 3 .and. has_latlon) then

         ilat = find_name(dnames, trim(lat_dim))
         ilon = find_name(dnames, trim(lon_dim))

         if (ilat > 0 .and. ilon > 0 .and. &
            ilat /= ilon .and. &
            itime /= ilat .and. itime /= ilon) then

            uses_latlon = .true.
            ok = .true.
            return
         end if
      end if

      errmsg = 'Variable '//trim(vname)//' has unsupported dimensions.'

   end subroutine resolve_timeseries_dims


   subroutine resolve_profile_dims(db, vname, time_dim, depth_dim, has_latlon, lat_dim, lon_dim, &
                                   uses_latlon, itime, idepth, ilat, ilon, ok, errmsg)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname, time_dim, depth_dim
      logical,      intent(in) :: has_latlon
      character(*), intent(in) :: lat_dim, lon_dim

      logical,      intent(out) :: uses_latlon
      integer,      intent(out) :: itime, idepth, ilat, ilon
      logical,      intent(out) :: ok
      character(*), intent(out) :: errmsg

      character(:), allocatable :: dnames(:)
      integer, allocatable :: dlens(:)
      integer :: ndims

      ok = .false.
      errmsg = ''

      uses_latlon = .false.
      itime  = 0
      idepth = 0
      ilat   = 0
      ilon   = 0

      call nc_var_dims(db, vname, dnames, dlens)
      ndims = size(dnames)

      itime  = find_name(dnames, trim(time_dim))
      idepth = find_name(dnames, trim(depth_dim))

      if (itime <= 0) then
         errmsg = 'Profile variable '//trim(vname)//' does not contain time dimension '// &
                  trim(time_dim)//'.'
         return
      end if

      if (idepth <= 0) then
         errmsg = 'Profile variable '//trim(vname)//' does not contain depth dimension '// &
                  trim(depth_dim)//'.'
         return
      end if

      if (itime == idepth) then
         errmsg = 'Profile variable '//trim(vname)//' uses the same dimension for time and depth.'
         return
      end if

      ! Point profile: var(time, depth), in either order.
      if (ndims == 2) then
         ok = .true.
         return
      end if

      ! Gridded profile: var(time, depth, lat, lon), in arbitrary order.
      if (ndims == 4 .and. has_latlon) then
         ilat = find_name(dnames, trim(lat_dim))
         ilon = find_name(dnames, trim(lon_dim))

         if (ilat > 0 .and. ilon > 0 .and. &
            ilat /= ilon .and. itime /= ilat .and. itime /= ilon .and. &
            idepth /= ilat .and. idepth /= ilon) then
            uses_latlon = .true.
            ok = .true.
            return
         end if
      end if

      errmsg = 'Profile variable '//trim(vname)//' has unsupported dimensions.'
   end subroutine resolve_profile_dims


   subroutine prepare_profile_spec(db, spec, time_dim, has_latlon, lat_dim, lon_dim, ok, errmsg)
      type(NcFile),   intent(in)    :: db
      type(DataSpec), intent(inout) :: spec
      character(*),   intent(in)    :: time_dim, lat_dim, lon_dim
      logical,        intent(in)    :: has_latlon
      logical,        intent(out)   :: ok
      character(*),   intent(out)   :: errmsg

      character(:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)
      logical :: uses_latlon, ok_dims
      integer :: itime, idepth, ilat, ilon, ndepth
      character(len=512) :: errmsg_dims
      real(rk), allocatable :: dz(:)

      ok = .false.
      errmsg = ''

      if (.not. allocated(spec%source_var)) then
         errmsg = 'Profile input is missing source_var.'
         return
      end if
      if (len_trim(spec%source_var) == 0) then
         errmsg = 'Profile input is missing source_var.'
         return
      end if

      if (.not. nc_has_var(db, trim(spec%source_var))) then
         errmsg = 'Profile variable '//trim(spec%source_var)//' is missing from NetCDF file.'
         return
      end if

      if (.not. allocated(spec%depth_var)) then
         errmsg = 'Profile variable '//trim(spec%source_var)//' requires depth_var.'
         return
      end if
      if (len_trim(spec%depth_var) == 0) then
         errmsg = 'Profile variable '//trim(spec%source_var)//' requires depth_var.'
         return
      end if

      if (.not. nc_has_var(db, trim(spec%depth_var))) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' for profile '// &
                  trim(spec%source_var)//' is missing from NetCDF file.'
         return
      end if

      call nc_var_dims(db, trim(spec%depth_var), dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' must be 1-D.'
         return
      end if

      spec%depth_dim = trim(dimnames(1))
      ndepth = dimlens(1)
      if (ndepth < 2) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' must contain at least two levels.'
         return
      end if

      if (trim(spec%depth_dim) == trim(time_dim)) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' cannot use the time dimension.'
         return
      end if

      if (has_latlon) then
         if (trim(spec%depth_dim) == trim(lat_dim) .or. trim(spec%depth_dim) == trim(lon_dim)) then
            errmsg = 'Depth coordinate '//trim(spec%depth_var)//' must use a dimension distinct from latitude/longitude.'
            return
         end if
      end if

      if (allocated(spec%source_depth)) deallocate(spec%source_depth)
      allocate(spec%source_depth(ndepth))
      call nc_read_real_1d(db, trim(spec%depth_var), spec%source_depth)

      if (any(.not. ieee_is_finite(spec%source_depth))) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' contains NaN/Inf values.'
         return
      end if

      allocate(dz(ndepth - 1))
      dz = spec%source_depth(2:ndepth) - spec%source_depth(1:ndepth - 1)
      if (.not. (all(dz > 0.0_rk) .or. all(dz < 0.0_rk))) then
         errmsg = 'Depth coordinate '//trim(spec%depth_var)//' must be strictly monotonic.'
         return
      end if

      call resolve_profile_dims(db, trim(spec%source_var), trim(time_dim), trim(spec%depth_dim), &
                                has_latlon, trim(lat_dim), trim(lon_dim), uses_latlon, &
                                itime, idepth, ilat, ilon, ok_dims, errmsg_dims)
      if (.not. ok_dims) then
         errmsg = trim(errmsg_dims)
         return
      end if

      ok = .true.
   end subroutine prepare_profile_spec


   subroutine check_var_dims_cf(db, vname, has_latlon, time_dim, lat_dim, lon_dim, present, missing)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname
      logical,      intent(in) :: has_latlon
      character(*), intent(in) :: time_dim, lat_dim, lon_dim
      character(:), allocatable, intent(inout) :: present(:), missing(:)

      logical :: uses_latlon, ok_dims
      integer :: itime, ilat, ilon
      character(len=512) :: errmsg_dims

      if (.not. nc_has_var(db, vname)) then
         call append_string(missing, vname)
         return
      end if

      call resolve_timeseries_dims(db, vname, time_dim, has_latlon, lat_dim, lon_dim, &
                                    uses_latlon, itime, ilat, ilon, ok_dims, errmsg_dims)

      if (.not. ok_dims) then
         call append_string(missing, vname)
         return
      end if

      call append_string(present, vname)

   end subroutine check_var_dims_cf


   logical function is_var_valid_in_period(db, vname, time_dim, lat_dim, lon_dim, has_latlon, &
                                           i0, i1, yi, xi) result(good)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname, time_dim, lat_dim, lon_dim
      logical,      intent(in) :: has_latlon
      integer,      intent(in) :: i0, i1, yi, xi

      real(rk), allocatable :: buf(:)
      integer :: nt

      good = .true.
      nt = max(0, i1 - i0 + 1)
      if (nt == 0) then
         good = .false.
         return
      end if

      allocate(buf(nt))
      call read_netcdf_timeseries_at_point(db, trim(vname), trim(time_dim), i0, i1, &
                                          has_latlon, trim(lat_dim), trim(lon_dim), &
                                          yi, xi, buf)
      if (any(.not. ieee_is_finite(buf))) good = .false.
   end function is_var_valid_in_period

   subroutine clear_scan(scan)
      type(NetcdfScan), intent(inout) :: scan

      if (allocated(scan%path)) deallocate(scan%path)
      if (allocated(scan%time_name)) deallocate(scan%time_name)
      if (allocated(scan%time_dim))  deallocate(scan%time_dim)
      if (allocated(scan%lat_name)) deallocate(scan%lat_name)
      if (allocated(scan%lat_dim))  deallocate(scan%lat_dim)
      if (allocated(scan%lon_name)) deallocate(scan%lon_name)
      if (allocated(scan%lon_dim))  deallocate(scan%lon_dim)
      if (allocated(scan%lat)) deallocate(scan%lat)
      if (allocated(scan%lon)) deallocate(scan%lon)
      if (allocated(scan%lon_convention)) deallocate(scan%lon_convention)
      if (allocated(scan%axis%t_s)) deallocate(scan%axis%t_s)
      if (allocated(scan%vars_present)) deallocate(scan%vars_present)
      if (allocated(scan%vars_missing)) deallocate(scan%vars_missing)

      scan%has_latlon = .false.
      scan%is_point = .true.
      scan%yi = 1
      scan%xi = 1
      scan%lat_found = 0.0_rk
      scan%lon_found = 0.0_rk
      scan%i0 = 1
      scan%i1 = 1
      scan%median_dt = 0.0_rk
      scan%is_regular = .true.
      scan%rel_max_dev = 0.0_rk
      scan%sim_offset = 0_lk
      scan%y0 = 0; scan%mon0 = 0; scan%d0 = 0; scan%h0 = 0; scan%mi0 = 0; scan%s0 = 0
      scan%y1 = 0; scan%mon1 = 0; scan%d1 = 0; scan%h1 = 0; scan%mi1 = 0; scan%s1 = 0
   end subroutine clear_scan


   subroutine clear_char_list(a)
      character(:), allocatable, intent(inout) :: a(:)
      if (allocated(a)) deallocate(a)
   end subroutine clear_char_list

end module data_loader_netcdf
