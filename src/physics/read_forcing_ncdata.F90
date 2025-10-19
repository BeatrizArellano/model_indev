module read_forcing_ncdata
  use precision_types, only: rk
  use netcdf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use netcdf_io,       only: NcFile, nc_open, nc_close, nc_has_var, nc_check, &
                             nc_var_dims, nc_read_real_1d, nc_get_att_str, nc_read_real_1d_slice

  use cf_time_utils,   only: CFUnits, TimeAxis, parse_cf_time,        &
                               seconds_since_datetime
  use calendar_types,  only: CFCalendar
  use time_utils,      only: detect_frequency, index_at_or_before, check_time_monotonic
  use str_utils,       only: inttostr, realtostr, list_to_str, append_string
  use find_utils,      only: find_name, has_name, argmin_abs_vec
  use geo_utils,       only: simple_distance_deg

  implicit none
  private

  public :: ForcingScan, scan_forcing, get_ts_slice_from_point

  !---------------------------------------------------------------------------
  ! Results/handle of Scan storing results from scanning a forcing NetCDF file.
  ! - Contains coordinate metadata (time, lat, lon) and validation results.
  ! - Record which requested variables are present/missing/derivable.
  ! - Time index range [i0:i1] covering the simulation time.
  !-------------------------------------------------------------------------
  type :: ForcingScan
     ! axis names
     character(:), allocatable :: time_name              ! Name used for the time coordinate
     character(:), allocatable :: lat_name, lon_name     ! e.g. 'lat', 'lon'

     ! geo
     logical :: has_latlon = .false.                     ! true if both lat & lon exist in file
     logical :: is_point   = .true.                      ! true if no lat/lon axes (t-only forcing)
     integer :: yi = 1, xi = 1                           ! nearest gridpoint indices to (lat0,lon0)
     real(rk), allocatable :: lat(:), lon(:)             ! 1-D coordinate vectors 
     character(:), allocatable :: lon_convention         ! '0_360' or '-180_180' or '-' if no lon
     real(rk) :: lat_found = 0.0_rk, lon_found = 0.0_rk  ! coords at (yi,xi), or requested if point

     ! time axis
     type(CFUnits)    :: u                               ! parsed time-units
     type(CFCalendar) :: cal                             ! parsed/derived calendar
     type(TimeAxis)   :: axis                            ! Contains t_s(:) [seconds since epoch], t_first, t_last
     integer :: i0 = 1, i1 = 1                           ! Time indices covering simulation time
     real(rk) :: step_s = 0.0_rk                         ! Time-step (median time-srep) in seconds
     logical  :: is_regular = .true.                     ! true if time-step deviations ≤ tolerance for all time-steps
     real(rk) :: rel_max_dev = 0.0_rk                    ! Maximum relative deviation (divided by time-step)

     ! variables
     character(:), allocatable :: vars_present(:)        ! Required variables found
     character(:), allocatable :: vars_missing(:)        ! required variables missing or with wrong dims
     character(:), allocatable :: vars_derivable(:)      ! Can be derived (e.g., wind_speed from uas,vas)

     ! Simulation time interval 
     integer :: y0, mon0, d0, h0, mi0, s0                ! start datetime (independent from reference date) UTC tieme
     integer :: y1, mon1, d1, h1, mi1, s1                ! end datetime for simulation
  end type ForcingScan

contains

   !=============================================================================
   ! Scans a forcing NetCDF file and validate that it contains the data
   ! for the site and the coverage for the simulation period. 
   ! Checks that the time coordinate is monotonic, and if lat/lon are present, 
   ! it finds the closest coordinate.
   ! - filename: path to NetCDF file
   ! - lat0, lon0: Site coordinates
   ! - y0..s0, y1..s1    [in]  Start and end UTC datetimes (year,mon,day,hour,min,sec).
   ! - calendar_default: When time:calendar missing (e.g., "gregorian")
   ! - required_vars: list of variable names to validate
   !   info:        ForcingScan handle with discovered metadata and indices.
   !   ok           true if scan succeeded; false otherwise.
   !   errmsg       Error message in case of any error
   !==============================================================================
   subroutine scan_forcing(filename, lat0, lon0,       &
                           y0,mon0,d0,h0,mi0,s0,       &
                           y1,mon1,d1,h1,mi1,s1,       &
                           calendar_default,required_vars, &
                           info, ok, errmsg, max_sep_deg)

      character(*),               intent(in)  :: filename
      real(rk),                   intent(in)  :: lat0, lon0
      integer,                    intent(in)  :: y0,mon0,d0,h0,mi0,s0 ! For initial simulation date
      integer,                    intent(in)  :: y1,mon1,d1,h1,mi1,s1 ! Date for end of simulation
      character(*),               intent(in)  :: calendar_default     ! Default calendar
      character(*), dimension(:), intent(in)  :: required_vars        ! List of required variables
      type(ForcingScan),          intent(out) :: info                 ! ForcingScan handle to request the information
      logical,                    intent(out) :: ok
      character(*),               intent(out) :: errmsg
      real(rk), optional,         intent(in)  :: max_sep_deg          ! Maximum allowed distance (deg) to grid-point

      type(NcFile) :: db
      character(len=:), allocatable :: time_name, latname, lonname
      character(len=:), allocatable :: cal_attr, units_attr
      logical :: has_time, has_lat, has_lon, pres
      logical :: has_ws, has_uas, has_vas, ws_required
      character(len=:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)
      integer :: ntime
      real(rk), allocatable :: time_orig(:)
      real(rk) :: t_start, t_end
      real(rk) :: lonq
      real(rk), parameter :: default_max_sep_deg = 0.25_rk
      real(rk) :: sep_limit, sep_deg
      logical  :: ok_parse
      integer :: i
      logical :: is_monotonic, has_equal_consecutive
      integer :: i_notmon

      ok = .false.; errmsg = ''

      ! Open netcdf file
      call nc_open(db, filename)

      ! Find time variable (expecting time for now)
      time_name = ''
      if (nc_has_var(db, 'time')) then
         time_name = 'time'; has_time = .true.
      else
         has_time = .false.
      end if
      if (.not. has_time) then
         errmsg = 'Missing required coordinate variable "time".'
         call nc_close(db); return
      end if
      info%time_name = time_name

      ! Get CF attributes (units required, calendar optional with default)
      call nc_get_att_str(db, time_name, 'units',   units_attr, pres)
      if (.not. pres) then
         errmsg = 'CF time: missing "units" attribute on variable '//trim(time_name)//'.'
         call nc_close(db); return
      end if
      call nc_get_att_str(db, time_name, 'calendar', cal_attr, pres)
      if (.not. pres) cal_attr = ''  ! will use default below

      ! Parse CF time units and calendar
      call parse_cf_time(units_attr, cal_attr, trim(calendar_default), info%u, info%cal, ok_parse, errmsg)
      if (.not. ok_parse) then
         call nc_close(db); return
      end if

      ! Read time dimension (1D, not spiral yet :))
      call nc_var_dims(db, time_name, dimnames, dimlens)
      if (size(dimlens) /= 1) then
         errmsg = 'CF time: variable '//trim(time_name)//' must be 1D.'; call nc_close(db); return
      end if
      ntime = dimlens(1)
      allocate(time_orig(ntime))
      call nc_read_real_1d(db, time_name, time_orig)

      ! Build TimeAxis using seconds since reference date
      info%axis%cal = info%cal
      info%axis%u   = info%u
      allocate(info%axis%t_s(ntime))
      info%axis%t_s    = time_orig * info%u%timeunit_to_seconds   !Secnds since reference datetime
      info%axis%t_first = info%axis%t_s(1)
      info%axis%t_last  = info%axis%t_s(ntime)


      ! Guard against NaN/Inf in time axis
      if (any(.not. ieee_is_finite(info%axis%t_s))) then
         errmsg = 'The time coordinate contains NaN/Inf values.'
         call nc_close(db); return
      end if
      ! Check whether the time-series is monotonic or has duplicate timestamps
      call check_time_monotonic(info%axis%t_s, is_monotonic, has_equal_consecutive, i_notmon)
      if (.not. is_monotonic) then
         errmsg = 'CF time is not sorted: First inversion at index ' // trim(adjustl(inttostr(i_notmon))) // &
                  '-' // trim(adjustl(inttostr(i_notmon+1))) // '). Please sort by time and reorder variables consistently.'
         call nc_close(db); return
      end if

      if (has_equal_consecutive) then
         errmsg = 'CF time contains duplicate timestamps. This is not supported.'
         call nc_close(db); return
      end if

      ! Detecting time-step in forcing data and regularity
      call detect_frequency(info%axis%t_s, info%step_s, info%is_regular, info%rel_max_dev)
      if (info%step_s <= 0.0_rk) then
         errmsg = 'Detected non-positive time step.'
         call nc_close(db); return
      end if


      !---- Find Latitude/Longitude --------------------------------- 
      !  - If present, read lat and lon, validate they’re not Nan or Inf.
      !  - Determine longitude convention ('0_360' vs '-180_180') 
      !  - Pick nearest gridpoint indices (yi, xi) and save the found coordinates.
      !  - If no lat/lon axes exist, we assume file as the forcing (is_point=true).
      latname = ''; lonname = ''
      has_lat = nc_has_var(db,'latitude') .or. nc_has_var(db,'lat')
      has_lon = nc_has_var(db,'longitude') .or. nc_has_var(db,'lon') .or. nc_has_var(db,'long')
      if (nc_has_var(db,'latitude')) then
         latname = 'latitude'                                              ! different names for latitude
      else
         latname = 'lat'
      end if
      if (has_lon) then
         if     (nc_has_var(db,'longitude')) then; lonname = 'longitude' ! Different names for longitude
         elseif (nc_has_var(db,'lon'))       then; lonname = 'lon'
         else                                    ; lonname = 'long'
         end if
      end if

      info%lat_name = latname            ! Assignas Name for the lat/lon coordinates
      info%lon_name = lonname

      info%has_latlon = (has_lat .and. has_lon); info%is_point = .not. info%has_latlon
      if (info%has_latlon) then
         ! --- 1-D vectors (no curvilinear support) ---
         call nc_var_dims(db, latname, dimnames, dimlens)
         if (size(dimlens) /= 1) then
            errmsg = 'Latitude must be 1D (curvilinear grids not supported).'
            call nc_close(db); return
         end if
         allocate(info%lat(dimlens(1)))
         call nc_read_real_1d(db, latname, info%lat)

         call nc_var_dims(db, lonname, dimnames, dimlens)
         if (size(dimlens) /= 1) then
            errmsg = 'Longitude must be 1D (curvilinear grids not supported).'
            call nc_close(db); return
         end if
         allocate(info%lon(dimlens(1)))
         call nc_read_real_1d(db, lonname, info%lon)

         ! --- Looking for NaN/Inf  ---
         if (any(.not. ieee_is_finite(info%lat)) .or. any(.not. ieee_is_finite(info%lon))) then
            errmsg = 'Latitude/Longitude contain NaN/Inf values.'
            call nc_close(db); return
         end if

         ! --- Determine longitude convention
         if (minval(info%lon) >= 0.0_rk .and. maxval(info%lon) <= 360.0_rk) then
            info%lon_convention = '0_360'
            lonq = modulo(lon0, 360.0_rk); if (lonq < 0.0_rk) lonq = lonq + 360.0_rk
         else
            info%lon_convention = '-180_180'
            lonq = modulo(lon0 + 180.0_rk, 360.0_rk) - 180.0_rk  ! (-180,180]
         end if

         ! Finding the nearest coordinates to the requested site
         info%yi = argmin_abs_vec(info%lat - lat0)
         info%xi = argmin_abs_vec(info%lon - lonq)
         info%lat_found = info%lat(info%yi)
         info%lon_found = info%lon(info%xi)

         if (present(max_sep_deg)) then
            sep_limit = max_sep_deg
         else
            sep_limit = default_max_sep_deg
         end if
         ! Check whether the site falls within max_sep_deg (default 0.25°)
         sep_deg = simple_distance_deg(lat0, lonq, info%lat_found, info%lon_found)
         if (sep_deg > sep_limit) then
            errmsg = 'No nearby grid point: nearest at ~'//realtostr(sep_deg,3)// &
                     '° (> '//realtostr(sep_limit,3)//'°) from requested site.'
            call nc_close(db); return
         end if

      else
         info%lon_convention = '-'
         info%lat_found = lat0; info%lon_found = lon0
         info%yi = 1; info%xi = 1
      end if


      !-----------Check coverage for simulation time---------------------------
      !Assign simulation datetimes
      info%y0=y0; info%mon0=mon0; info%d0=d0; info%h0=h0; info%mi0=mi0; info%s0=s0
      info%y1=y1; info%mon1=mon1; info%d1=d1; info%h1=h1; info%mi1=mi1; info%s1=s1

      ! Convert requested start/end datetimes to seconds since reference date
      t_start = seconds_since_datetime(info%cal, info%u, y0, mon0, d0, h0, mi0, s0)
      t_end   = seconds_since_datetime(info%cal, info%u, y1, mon1, d1, h1, mi1, s1)
      ! Coverage check
      if (t_start < info%axis%t_first .or. t_end > info%axis%t_last .or. t_end < t_start) then
         errmsg = 'Forcing data does not cover simulation period'
         call nc_close(db); return
      end if

      ! Indices for [t_start, t_end]
      info%i0 = index_at_or_before(info%axis%t_s, t_start)
      info%i1 = index_at_or_before(info%axis%t_s, t_end)
      if (info%i0 < 1) info%i0 = 1
      if (info%i1 < info%i0) info%i1 = info%i0

      !---- Variable bookkeeping -----------------------------------------------------
      ! Validate required variables and that they have required coordinates.
      ! A variable is as existing if it exists and its dimensions include:
      !  time or lat/lon (if present)
      call clear_char_lists(info%vars_present, info%vars_missing, info%vars_derivable)
      do i=1, size(required_vars)
         call check_var_dims_cf(db, trim(required_vars(i)), info%has_latlon, info%time_name, info%lat_name, info%lon_name, &
                                 info%vars_present, info%vars_missing)
      end do
      ! Derivable: wind_speed from uas, vas (only if requested)
      has_ws = allocated(info%vars_present)  .and. has_name(info%vars_present,'wind_speed')
      has_uas = allocated(info%vars_present) .and. has_name(info%vars_present,'uas')
      has_vas = allocated(info%vars_present) .and. has_name(info%vars_present,'vas')
      ws_required = has_name(required_vars,'wind_speed')
      if (.not. has_ws) then
         if (has_uas .and. has_vas) then
            call append_string(info%vars_derivable, 'wind_speed')
         elseif (ws_required) then
            call append_string(info%vars_missing, 'wind_speed')
         end if
      end if
      ! If any required variable is missing fail with message.
      if (allocated(info%vars_missing) .and. size(info%vars_missing) > 0) then
         errmsg = 'Missing required variables: ' // trim(list_to_str(info%vars_missing, ', '))
         call nc_close(db); return
      end if

      !---- Final data-quality check: ensure Variables have valid values during simulation period ----
      if (allocated(info%vars_present)) then
         do i = 1, size(info%vars_present)
            if (.not. is_var_valid_in_period(db, trim(info%vars_present(i)), info%time_name, &
                                             info%lat_name, info%lon_name, info%has_latlon, &
                                             info%i0, info%i1, info%yi, info%xi)) then
               errmsg = 'Variable "'//trim(info%vars_present(i))// &
                        '" contains NaN/Inf in the requested simulation period at the selected site.'
               call nc_close(db); return
            end if
         end do
      end if

      ! All checks passed and indices computed. The ForcingScan handle 'info' is ready
      ok = .true.
      call nc_close(db)   ! Close file
   
   end subroutine scan_forcing


   !------------------------------------------------------------------------------
   ! Returns a slice of a time series between the indices i0 and i1
   ! at a grid-point (yi,xi) corresponding to the requested variable
   ! that has dimensions (time) or (time,lat,lon) in ANY order.
   !
   ! Args:
   !   db          [in]  NcFile handle (already open)
   !   varname     [in]  Variable to read (e.g., 'shortwave')
   !   time_name   [in]  Time dimension name discovered by scan ('time')
   !   i0,i1       [in]  Indices along time (i1 >= i0)
   !   has_latlon  [in]  .true. if the var is on a lat/lon grid (time,lat,lon)
   !   lat_name    [in]  Latitude dimension name (ignored if has_latlon=.false.)
   !   lon_name    [in]  Longitude dimension name (ignored if has_latlon=.false.)
   !   yi, xi      [in]  Selected gridpoint indices (ignored if has_latlon=.false.)
   !   out(:)      [out] Output buffer sized to i1-i0+1 (real(rk))
   !
   !------------------------------------------------------------------------------
   subroutine get_ts_slice_from_point(db, varname, time_name, i0, i1, &
                                      has_latlon, lat_name, lon_name, yi, xi, out)
  
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: varname, time_name, lat_name, lon_name
      logical,      intent(in) :: has_latlon
      integer,      intent(in) :: i0, i1, yi, xi
      real(rk),     intent(out):: out(:)

      integer :: vid, ndims, dimids(NF90_MAX_VAR_DIMS), xtype, natts
      character(len=:), allocatable :: dnames(:)
      integer, allocatable :: dlens(:)
      integer :: itime, ilat, ilon, nt
      integer :: start(NF90_MAX_VAR_DIMS), count(NF90_MAX_VAR_DIMS)

      call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid),        'inq_varid('//trim(varname)//')')
      call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                     'inquire_variable('//trim(varname)//')')
      if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'timeseries_at_point: '//trim(varname)//' is character')

      call nc_var_dims(db, varname, dnames, dlens)

      itime = find_name(dnames, trim(time_name))
      if (itime <= 0) call nc_check(NF90_EBADDIM, 'timeseries_at_point: time dim not found for '//trim(varname))

      if (has_latlon) then
         ilat = find_name(dnames, trim(lat_name))
         ilon = find_name(dnames, trim(lon_name))
         if (ilat <= 0 .or. ilon <= 0) &
            call nc_check(NF90_EBADDIM, 'timeseries_at_point: lat/lon dims not found for '//trim(varname))
         if (yi < 1 .or. yi > dlens(ilat) .or. xi < 1 .or. xi > dlens(ilon)) &
            call nc_check(NF90_EEDGE, 'timeseries_at_point: yi/xi out of range for '//trim(varname))
      else
         if (ndims /= 1) call nc_check(NF90_EINVALCOORDS, 'timeseries_at_point: '//trim(varname)//' must be 1-D (time)')
      end if

      nt = i1 - i0 + 1
      if (nt < 1)                           call nc_check(NF90_EEDGE, 'timeseries_at_point: empty time window for '//trim(varname))
      if (i0 < 1 .or. i1 > dlens(itime))    call nc_check(NF90_EEDGE, 'timeseries_at_point: time indices out of range for '//trim(varname))
      if (size(out) /= nt)                  call nc_check(NF90_EEDGE, 'timeseries_at_point: output size mismatch for '//trim(varname))

      start(1:ndims) = 1;  count(1:ndims) = 1
      start(itime)   = i0; count(itime)   = nt
      if (has_latlon) then
         start(ilat) = yi; count(ilat)    = 1
         start(ilon) = xi; count(ilon)    = 1
      end if

      call nc_check(nf90_get_var(db%ncid, vid, out, start=start(1:ndims), count=count(1:ndims)), &
                     'get_var slice '//trim(varname))
   end subroutine get_ts_slice_from_point

  !------------------------
  ! Local Helpers
  !------------------------


   !------------------------------------------------------------------------------
   ! Validates that variable vname exists and has the expected coordinates
   !   - Requires the time dimension named time_name.
   !   - If require_latlon=true, also requires lat_name and lon_name.
   !
   !   - If the variable is missing OR any required axis name is absent from its
   !     dimension list, the variable name is declared as missing.
   !   - Otherwise, it is declared as present.
   subroutine check_var_dims_cf(db, vname, require_latlon, time_name, lat_name, lon_name, present, missing)
      type(NcFile), intent(in) :: db
      character(*), intent(in) :: vname
      logical, intent(in) :: require_latlon
      character(*), intent(in) :: time_name
      character(*), intent(in) :: lat_name, lon_name
      character(:), allocatable, intent(inout) :: present(:), missing(:)
      character(len=:), allocatable :: dimnames(:)
      integer, allocatable :: dimlens(:)

      if (.not. nc_has_var(db, vname)) then
         call append_string(missing, vname); return
      end if

      call nc_var_dims(db, vname, dimnames, dimlens)
      if (.not. has_name(dimnames, trim(time_name))) then
         call append_string(missing, vname); return
      end if
      if (require_latlon) then
         if (.not. has_name(dimnames, trim(lat_name))) then
            call append_string(missing, vname); return
         end if
         if (.not. has_name(dimnames, trim(lon_name))) then
            call append_string(missing, vname); return
         end if
      end if
      call append_string(present, vname)
   end subroutine check_var_dims_cf

   ! Deallocates up to three allocatable character arrays
   subroutine clear_char_lists(a,b,c)
      character(:), allocatable, intent(inout), optional :: a(:), b(:), c(:)
      if (present(a)) then; if (allocated(a)) deallocate(a); end if
      if (present(b)) then; if (allocated(b)) deallocate(b); end if
      if (present(c)) then; if (allocated(c)) deallocate(c); end if
   end subroutine clear_char_lists 
   
   
   ! Verifies that the requested variable contains valid numeric data (Not Nan or Inf) within a time-period
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
         good = .false.; return
      end if
      allocate(buf(nt))
      call get_ts_slice_from_point(db, vname, time_name, i0, i1, has_latlon, lat_name, lon_name, yi, xi, buf)
      if (any(.not. ieee_is_finite(buf))) good = .false.
      deallocate(buf)
   end function is_var_valid_in_period
  
end module read_forcing_ncdata
