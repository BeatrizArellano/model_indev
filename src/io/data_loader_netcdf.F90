module data_loader_netcdf
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_types, only: rk, lk
   use netcdf
   use netcdf_io, only: NcFile, nc_open, nc_close, nc_has_var, nc_check, &
                        nc_var_dims, nc_read_real_1d, nc_get_att_str
   use data_types, only: DataSpec, DataVarSeries, DATA_INPUT_FILE
   use time_types, only: CFUnits, TimeAxis, CFCalendar, DateTime
   use cf_time_utils, only: parse_cf_time, seconds_since_datetime_file
   use time_utils, only: detect_frequency, index_at_or_before, index_at_or_after, check_time_monotonic
   use str_utils, only: inttostr, realtostr, list_to_str, append_string
   use find_utils, only: find_name, has_name, argmin_abs_vec
   use geo_utils, only: LocationInfo, simple_distance_deg

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
      character(:), allocatable :: lat_name
      character(:), allocatable :: lon_name

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
      character(*),        intent(in)  :: path
      type(DataSpec),      intent(in)  :: specs(:)
      type(LocationInfo),  intent(in)  :: location
      type(DateTime),      intent(in)  :: start_datetime, end_datetime
      character(*),        intent(in)  :: calendar_default
      type(NetcdfScan),    intent(out) :: scan
      logical,             intent(out) :: ok
      character(*),        intent(out) :: errmsg
      real(rk), optional,  intent(in)  :: max_sep_deg

      type(NcFile) :: db
      character(:), allocatable :: required_vars(:)
      character(:), allocatable :: time_name
      character(:), allocatable :: lat_name, lon_name
      character(:), allocatable :: units_attr, cal_attr
      character(:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)
      real(rk), allocatable :: time_orig(:)
      logical :: pres, ok_parse
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

      call find_horizontal_coordinates(db, location, scan, ok, errmsg, max_sep_deg)
      if (.not. ok) then
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

      scan%sim_offset = int(nint(t_start), lk)
      scan%i0 = max(1, index_at_or_before(scan%axis%t_s, t_start))
      scan%i1 = max(scan%i0, index_at_or_before(scan%axis%t_s, t_end))

      call clear_char_list(scan%vars_present)
      call clear_char_list(scan%vars_missing)

      do i = 1, size(required_vars)
         call check_var_dims_cf(db, trim(required_vars(i)), scan%has_latlon, scan%time_name, &
                                scan%lat_name, scan%lon_name, scan%vars_present, scan%vars_missing)
      end do

      if (allocated(scan%vars_missing) .and. size(scan%vars_missing) > 0) then
         errmsg = 'Missing or incompatible NetCDF variables in '//trim(path)//': '// &
                  trim(list_to_str(scan%vars_missing, ', '))
         call nc_close(db)
         return
      end if

      do i = 1, size(required_vars)
         if (.not. is_var_valid_in_period(db, trim(required_vars(i)), scan%time_name, &
                                          scan%lat_name, scan%lon_name, scan%has_latlon, &
                                          scan%i0, scan%i1, scan%yi, scan%xi)) then
            errmsg = 'Variable '//trim(required_vars(i))// &
                     ' contains NaN/Inf values over the requested interval in '//trim(path)//'.'
            call nc_close(db)
            return
         end if
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

      series%t_axis = int(nint(scan%axis%t_s(i0:i1)), lk)

      call read_netcdf_timeseries_at_point(db, trim(spec%source_var), trim(scan%time_name), i0, i1, &
                                           scan%has_latlon, trim(scan%lat_name), trim(scan%lon_name), &
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
         dt_last = max(1_lk, int(nint(scan%median_dt), lk))

         series%t_edge(1) = series%t_axis(1)
         series%t_edge(2) = series%t_axis(1) + dt_last

         series%t_next = INF_EDGE
      end if
   end subroutine load_netcdf_series


   subroutine read_netcdf_timeseries_at_point(db, varname, time_name, i0, i1, &
                                              has_latlon, lat_name, lon_name, yi, xi, out)
      type(NcFile), intent(in)  :: db
      character(*), intent(in)  :: varname, time_name, lat_name, lon_name
      logical,      intent(in)  :: has_latlon
      integer,      intent(in)  :: i0, i1, yi, xi
      real(rk),     intent(out) :: out(:)

      integer :: vid, ndims, dimids(NF90_MAX_VAR_DIMS), xtype, natts
      character(:), allocatable :: dnames(:)
      integer, allocatable :: dlens(:)
      integer :: itime, ilat, ilon, nt
      integer :: start(NF90_MAX_VAR_DIMS), count(NF90_MAX_VAR_DIMS)

      call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
      call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                    'inquire_variable('//trim(varname)//')')

      if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'read_netcdf_timeseries_at_point: '//trim(varname)//' is character')

      call nc_var_dims(db, varname, dnames, dlens)

      itime = find_name(dnames, trim(time_name))
      if (itime <= 0) call nc_check(NF90_EBADDIM, 'time dimension not found for '//trim(varname))

      if (has_latlon) then
         ilat = find_name(dnames, trim(lat_name))
         ilon = find_name(dnames, trim(lon_name))
         if (ilat <= 0 .or. ilon <= 0) call nc_check(NF90_EBADDIM, 'lat/lon dimensions not found for '//trim(varname))
         if (yi < 1 .or. yi > dlens(ilat) .or. xi < 1 .or. xi > dlens(ilon)) &
            call nc_check(NF90_EEDGE, 'selected lat/lon index out of range for '//trim(varname))
      else
         if (ndims /= 1) call nc_check(NF90_EINVALCOORDS, trim(varname)//' must be 1-D when no lat/lon axes are present')
      end if

      nt = i1 - i0 + 1
      if (nt < 1) call nc_check(NF90_EEDGE, 'empty time window for '//trim(varname))
      if (i0 < 1 .or. i1 > dlens(itime)) call nc_check(NF90_EEDGE, 'time indices out of range for '//trim(varname))
      if (size(out) /= nt) call nc_check(NF90_EEDGE, 'output size mismatch for '//trim(varname))

      start(1:ndims) = 1
      count(1:ndims) = 1
      start(itime) = i0
      count(itime) = nt

      if (has_latlon) then
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
         if (.not. allocated(specs(i)%path)) cycle
         if (trim(specs(i)%path) /= trim(path)) cycle
         if (.not. allocated(specs(i)%source_var)) cycle
         maxlen = max(maxlen, len_trim(specs(i)%source_var))
      end do

      allocate(character(len=maxlen) :: names(size(specs)))
      names = ''

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
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
      allocate(scan%lat(dimlens(1)))
      call nc_read_real_1d(db, scan%lat_name, scan%lat)

      call nc_var_dims(db, scan%lon_name, dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'Longitude coordinate must be 1-D; curvilinear grids are not supported yet.'
         return
      end if
      allocate(scan%lon(dimlens(1)))
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


   subroutine check_var_dims_cf(db, vname, require_latlon, time_name, lat_name, lon_name, present, missing)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname
      logical,      intent(in) :: require_latlon
      character(*), intent(in) :: time_name, lat_name, lon_name
      character(:), allocatable, intent(inout) :: present(:), missing(:)

      character(:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)

      if (.not. nc_has_var(db, vname)) then
         call append_string(missing, vname)
         return
      end if

      call nc_var_dims(db, vname, dimnames, dimlens)

      if (.not. has_name(dimnames, trim(time_name))) then
         call append_string(missing, vname)
         return
      end if

      if (require_latlon) then
         if (.not. has_name(dimnames, trim(lat_name))) then
            call append_string(missing, vname)
            return
         end if

         if (.not. has_name(dimnames, trim(lon_name))) then
            call append_string(missing, vname)
            return
         end if
      end if

      call append_string(present, vname)
   end subroutine check_var_dims_cf


   logical function is_var_valid_in_period(db, vname, time_name, lat_name, lon_name, has_latlon, &
                                           i0, i1, yi, xi) result(good)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname, time_name, lat_name, lon_name
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
      call read_netcdf_timeseries_at_point(db, vname, time_name, i0, i1, has_latlon, lat_name, lon_name, yi, xi, buf)
      if (any(.not. ieee_is_finite(buf))) good = .false.
   end function is_var_valid_in_period


   subroutine clear_scan(scan)
      type(NetcdfScan), intent(inout) :: scan

      if (allocated(scan%path)) deallocate(scan%path)
      if (allocated(scan%time_name)) deallocate(scan%time_name)
      if (allocated(scan%lat_name)) deallocate(scan%lat_name)
      if (allocated(scan%lon_name)) deallocate(scan%lon_name)
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