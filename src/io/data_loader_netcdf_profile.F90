module data_loader_netcdf_profile
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   use data_loader_netcdf, only: NetcdfScan, read_netcdf_profile_at_point
   use data_types,         only: DataSpec, DataVarSeries
   use netcdf_io,          only: NcFile
   use precision_types,    only: rk, lk

   implicit none
   private

   public :: load_netcdf_profile_series

   integer(lk), parameter :: INF_EDGE = huge(1_lk)

contains

   subroutine load_netcdf_profile_series(series, scan, db, spec, i0, i1)
      type(DataVarSeries), intent(inout) :: series
      type(NetcdfScan),    intent(in)    :: scan
      type(NcFile),        intent(in)    :: db
      type(DataSpec),      intent(in)    :: spec
      integer,             intent(in)    :: i0, i1

      real(rk), allocatable :: native_values(:,:)
      integer :: nt, ntarget, it, i
      integer(lk) :: dt_last
      real(rk) :: source_min, source_max, target_min, target_max

      if (.not. spec%is_profile) then
         call profile_load_error('Spec is not marked as a profile for ', spec%name)
      end if

      if (.not. allocated(spec%source_depth)) then
         call profile_load_error('source_depth is not available for ', spec%name)
      end if
      if (size(spec%source_depth) < 2) then
         call profile_load_error('source_depth must contain at least two levels for ', spec%name)
      end if

      if (.not. allocated(spec%target_depth)) then
         call profile_load_error('target_depth is not available for ', spec%name)
      end if
      if (size(spec%target_depth) < 1) then
         call profile_load_error('target_depth is empty for ', spec%name)
      end if
      if (any(.not. ieee_is_finite(spec%target_depth))) then
         call profile_load_error('target_depth contains NaN/Inf for ', spec%name)
      end if

      nt = max(0, i1 - i0 + 1)
      if (nt <= 0) then
         call profile_load_error('Empty time window for ', spec%name)
      end if

      call read_netcdf_profile_at_point(db, trim(spec%source_var), trim(scan%time_dim), &
                                        trim(spec%depth_dim), i0, i1, scan%has_latlon, &
                                        trim(scan%lat_dim), trim(scan%lon_dim), &
                                        scan%yi, scan%xi, native_values)

      if (size(native_values, 1) /= size(spec%source_depth)) then
         call profile_load_error('Source depth size does not match profile data for ', spec%name)
      end if
      if (size(native_values, 2) /= nt) then
         call profile_load_error('Time size does not match profile data for ', spec%name)
      end if

      source_min = minval(spec%source_depth)
      source_max = maxval(spec%source_depth)
      target_min = minval(spec%target_depth)
      target_max = maxval(spec%target_depth)

      if (target_max < source_min .or. target_min > source_max) then
         call profile_load_error('Source and target depth ranges do not overlap for ', spec%name)
      end if

      if (allocated(series%t_axis)) deallocate(series%t_axis)
      if (allocated(series%t_edge)) deallocate(series%t_edge)
      if (allocated(series%values)) deallocate(series%values)
      if (allocated(series%depth)) deallocate(series%depth)
      if (allocated(series%profile_values)) deallocate(series%profile_values)

      ntarget = size(spec%target_depth)

      allocate(series%t_axis(nt))
      allocate(series%depth(ntarget))
      allocate(series%profile_values(ntarget, nt))

      series%name        = trim(spec%name)
      series%units       = trim(spec%units)
      series%is_const    = .false.
      series%is_profile  = .true.
      series%n           = nt
      series%idx         = 1
      series%depth       = spec%target_depth
      series%t_axis      = nint(scan%axis%t_s(i0:i1), kind=lk)

      do it = 1, nt
         call interpolate_profile_linear(spec%source_depth, native_values(:,it), &
                                         spec%target_depth, series%profile_values(:,it))
      end do

      allocate(series%t_edge(nt + 1))

      ! Same left/step-ahead time convention used by scalar NetCDF series:
      ! profile(:,i) is valid over [t_axis(i), t_axis(i+1)).
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
   end subroutine load_netcdf_profile_series


   subroutine interpolate_profile_linear(source_depth, source_values, target_depth, target_values)
      real(rk), intent(in)  :: source_depth(:)
      real(rk), intent(in)  :: source_values(:)
      real(rk), intent(in)  :: target_depth(:)
      real(rk), intent(out) :: target_values(:)

      integer :: nsource, ntarget, j, k
      logical :: increasing
      real(rk) :: z, w

      nsource = size(source_depth)
      ntarget = size(target_depth)

      if (nsource < 2) then
         error stop 'interpolate_profile_linear: source profile must contain at least two levels.'
      end if
      if (size(source_values) /= nsource) then
         error stop 'interpolate_profile_linear: source depth/value size mismatch.'
      end if
      if (size(target_values) /= ntarget) then
         error stop 'interpolate_profile_linear: target depth/value size mismatch.'
      end if

      increasing = source_depth(nsource) > source_depth(1)

      do k = 1, ntarget
         z = target_depth(k)

         if (increasing) then
            ! Do not linearly extrapolate beyond the source profile. Hold the
            ! nearest endpoint value instead.
            if (z <= source_depth(1)) then
               target_values(k) = source_values(1)
               cycle
            else if (z >= source_depth(nsource)) then
               target_values(k) = source_values(nsource)
               cycle
            end if

            do j = 1, nsource - 1
               if (z >= source_depth(j) .and. z <= source_depth(j + 1)) then
                  w = (z - source_depth(j)) / (source_depth(j + 1) - source_depth(j))
                  target_values(k) = source_values(j) + w * (source_values(j + 1) - source_values(j))
                  exit
               end if
            end do
         else
            if (z >= source_depth(1)) then
               target_values(k) = source_values(1)
               cycle
            else if (z <= source_depth(nsource)) then
               target_values(k) = source_values(nsource)
               cycle
            end if

            do j = 1, nsource - 1
               if (z <= source_depth(j) .and. z >= source_depth(j + 1)) then
                  w = (z - source_depth(j)) / (source_depth(j + 1) - source_depth(j))
                  target_values(k) = source_values(j) + w * (source_values(j + 1) - source_values(j))
                  exit
               end if
            end do
         end if
      end do
   end subroutine interpolate_profile_linear


   subroutine profile_load_error(message, name)
      character(*), intent(in) :: message, name

      write(*,'(A)') 'ERROR load_netcdf_profile_series: '//trim(message)//trim(name)
      error stop 'load_netcdf_profile_series failed.'
   end subroutine profile_load_error

end module data_loader_netcdf_profile
