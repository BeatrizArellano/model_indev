module data_loader
   use data_types,          only: DataLoaderCfg, DataSpec, DataFileInfo, DataVarSeries, InputData, &
                                  DATA_INPUT_FILE, DATA_INPUT_CONSTANT, DATA_INPUT_COMPUTE, DATA_INPUT_OFF, &
                                  DATA_FORMAT_UNKNOWN, DATA_FORMAT_NETCDF, DATA_FORMAT_CSV, &
                                  DATA_TIME_ABSOLUTE, DATA_TIME_REPEAT_YEAR
   use data_loader_netcdf,  only: NetcdfScan, scan_netcdf_file, &
                                  build_netcdf_year_windows, load_netcdf_series
   use geo_utils,           only: LocationInfo
   use netcdf_io,           only: NcFile, nc_open, nc_close
   use precision_types,     only: rk, lk
   use str_utils,           only: to_lower
   use time_types,          only: DateTime, CFCalendar, calendar_compatible, cal_unknown

   implicit none
   private

   public :: DataLoaderState
   public :: scan_and_init_data
   public :: load_input_data
   public :: print_data_summary
   public :: find_data_index
   public :: init_series_cursor
   public :: value_at_step

   integer(lk), parameter :: INF_EDGE = huge(1_lk)

   type :: DataLoaderState
      type(DataLoaderCfg) :: cfg

      type(DataSpec),     allocatable :: specs(:)
      type(DataFileInfo), allocatable :: files(:)
      type(NetcdfScan),   allocatable :: nc_scans(:)

      type(CFCalendar) :: sim_cal

      integer :: sim_y_start = 0
      integer :: sim_y_end   = 0
   end type DataLoaderState

contains

   subroutine scan_and_init_data(specs_in, cfg, calendar_cfg, location, start_datetime, end_datetime, state, ok, errmsg)
      type(DataSpec),      intent(in)  :: specs_in(:)
      type(DataLoaderCfg), intent(in)  :: cfg
      type(CFCalendar),    intent(in)  :: calendar_cfg
      type(LocationInfo),  intent(in)  :: location
      type(DateTime),      intent(in)  :: start_datetime, end_datetime
      type(DataLoaderState), intent(out) :: state
      logical,             intent(out) :: ok
      character(*),        intent(out) :: errmsg

      type(DataSpec), allocatable :: specs(:)
      type(DataFileInfo), allocatable :: files(:)

      type(DateTime) :: scan_start_datetime, scan_end_datetime

      integer :: i, j, ref_calendar
      logical :: lok
      character(len=512) :: lmsg

      ok = .false.
      errmsg = ''

      call clear_state(state)

      allocate(specs(size(specs_in)))
      specs = specs_in

      call normalise_specs(specs)
      call validate_specs(specs, ok, errmsg)
      if (.not. ok) return

      call build_file_list(specs, files)

      do i = 1, size(specs)
         if (specs(i)%input_type == DATA_INPUT_FILE) then
            specs(i)%file_index = find_file_index(files, specs(i)%path)

            if (specs(i)%file_index <= 0) then
               errmsg = 'scan_and_init_data: file not found for variable '//trim(specs(i)%name)
               return
            end if
         end if
      end do

      allocate(state%nc_scans(size(files)))

      do j = 1, size(files)
         select case (files(j)%format)

         case (DATA_FORMAT_NETCDF)
            if (cfg%repeat_enabled) then
               scan_start_datetime = DateTime(cfg%repeat_year, 1, 1, 0, 0, 0)
               scan_end_datetime   = DateTime(cfg%repeat_year + 1, 1, 1, 0, 0, 0)
            else
               scan_start_datetime = start_datetime
               scan_end_datetime   = end_datetime
            end if

            call scan_netcdf_file(files(j)%path, specs, location, scan_start_datetime, scan_end_datetime, &
                                 calendar_cfg%name(), state%nc_scans(j), lok, lmsg)

            if (.not. lok) then
               errmsg = trim(lmsg)
               return
            end if

         case (DATA_FORMAT_CSV)
            errmsg = 'CSV input files are recognised but not implemented yet: '//trim(files(j)%path)
            return

         case default
            errmsg = 'Unsupported input file format for '//trim(files(j)%path)
            return
         end select
      end do

      do i = 1, size(specs)
         if (specs(i)%input_type == DATA_INPUT_FILE) then
            j = specs(i)%file_index

            select case (files(j)%format)
            case (DATA_FORMAT_NETCDF)
               specs(i)%calendar   = state%nc_scans(j)%cal%kind
               specs(i)%first_time = state%nc_scans(j)%axis%t_first
               specs(i)%last_time  = state%nc_scans(j)%axis%t_last
            end select
         end if
      end do

      state%sim_y_start = start_datetime%year
      state%sim_y_end   = end_datetime%year

      ref_calendar = cal_unknown

      do j = 1, size(files)
         if (files(j)%format == DATA_FORMAT_NETCDF) then
            if (ref_calendar == cal_unknown) then
               ref_calendar = state%nc_scans(j)%cal%kind
            else if (.not. calendar_compatible(state%nc_scans(j)%cal%kind, ref_calendar)) then
               errmsg = 'Calendar mismatch across input files. Preprocess inputs to a common calendar.'
               return
            end if
         end if
      end do

      if (ref_calendar /= cal_unknown) then
         if (cfg%cfg_calendar == cal_unknown) then
            state%sim_cal%kind = ref_calendar
         else
            if (.not. calendar_compatible(cfg%cfg_calendar, ref_calendar)) then
               errmsg = 'Calendar mismatch between configuration and input data.'
               return
            end if

            state%sim_cal%kind = cfg%cfg_calendar
         end if
      else
         state%sim_cal%kind = cfg%cfg_calendar
      end if

      do i = 1, size(specs)
         if (specs(i)%input_type == DATA_INPUT_FILE) then
            j = specs(i)%file_index

            select case (files(j)%format)
            case (DATA_FORMAT_NETCDF)
               if (.not. cfg%load_yearly) then
                  call build_full_window(specs(i), state%nc_scans(j))
               else if (cfg%repeat_enabled) then
                  call build_netcdf_year_windows(specs(i), state%nc_scans(j), &
                                                 cfg%repeat_year, cfg%repeat_year, &
                                                 scan_start_datetime, scan_end_datetime)
               else
                  call build_netcdf_year_windows(specs(i), state%nc_scans(j), &
                                                 state%sim_y_start, state%sim_y_end, &
                                                 start_datetime, end_datetime)
               end if
            end select
         end if
      end do

      state%cfg = cfg
      if (state%cfg%repeat_enabled) state%cfg%time_mode = DATA_TIME_REPEAT_YEAR

      allocate(state%specs(size(specs)))
      state%specs = specs

      allocate(state%files(size(files)))
      state%files = files

      ok = .true.
   end subroutine scan_and_init_data


   subroutine load_input_data(state, year_k, input, ok, errmsg)
      type(DataLoaderState), intent(in)    :: state
      integer,               intent(in)    :: year_k
      type(InputData),       intent(inout) :: input
      logical,               intent(out)   :: ok
      character(*),          intent(out)   :: errmsg

      type(NcFile) :: db
      integer :: i, j, i0, i1, nt, win_k

      ok = .false.
      errmsg = ''

      call ensure_input_data_size(input, size(state%specs))

      do i = 1, size(state%specs)
         input%vars(i)%name = state%specs(i)%name
      end do

      do j = 1, size(state%files)
         select case (state%files(j)%format)

         case (DATA_FORMAT_NETCDF)
            call nc_open(db, state%files(j)%path)

            do i = 1, size(state%specs)
               if (state%specs(i)%input_type /= DATA_INPUT_FILE) cycle
               if (state%specs(i)%file_index /= j) cycle

               if (.not. allocated(state%specs(i)%idx_window)) then
                  errmsg = 'load_input_data: idx_window not built for '//trim(state%specs(i)%name)
                  call nc_close(db)
                  return
               end if

               if ((.not. state%cfg%load_yearly) .or. state%cfg%repeat_enabled) then
                  win_k = 1
               else
                  win_k = year_k
               end if

               if (win_k < 1 .or. win_k > size(state%specs(i)%idx_window, 2)) then
                  errmsg = 'load_input_data: year window index out of range for '//trim(state%specs(i)%name)
                  call nc_close(db)
                  return
               end if

               i0 = state%specs(i)%idx_window(1, win_k)
               i1 = state%specs(i)%idx_window(2, win_k)
               nt = max(0, i1 - i0 + 1)

               if (nt <= 0) then
                  errmsg = 'load_input_data: empty time slice for '//trim(state%specs(i)%name)
                  call nc_close(db)
                  return
               end if

               call load_netcdf_series(input%vars(i), state%nc_scans(j), db, state%specs(i), i0, i1)

               input%vars(i)%name  = state%specs(i)%name
               input%vars(i)%units = state%specs(i)%units

               input%vars(i)%cal         = state%nc_scans(j)%cal
               input%vars(i)%u           = state%nc_scans(j)%u
               input%vars(i)%sim_offset  = state%nc_scans(j)%sim_offset
               input%vars(i)%time_mode   = state%cfg%time_mode
               input%vars(i)%repeat_year = state%cfg%repeat_year
            end do

            call nc_close(db)

         case (DATA_FORMAT_CSV)
            errmsg = 'CSV loading is recognised but not implemented yet: '//trim(state%files(j)%path)
            return

         case default
            errmsg = 'Unsupported input file format for '//trim(state%files(j)%path)
            return
         end select
      end do

      do i = 1, size(state%specs)
         select case (state%specs(i)%input_type)

         case (DATA_INPUT_CONSTANT)
            call set_constant_series(input%vars(i), state%specs(i)%name, state%specs(i)%const_value)

         case (DATA_INPUT_OFF)
            call set_constant_series(input%vars(i), state%specs(i)%name, 0.0_rk)

         case (DATA_INPUT_FILE)
            ! Already loaded above.

         case (DATA_INPUT_COMPUTE)
            errmsg = 'load_input_data: DATA_INPUT_COMPUTE is not handled by the generic loader for '// &
                     trim(state%specs(i)%name)
            return

         case default
            errmsg = 'load_input_data: unsupported input_type for '//trim(state%specs(i)%name)
            return
         end select
      end do

      do i = 1, size(input%vars)
         call init_series_cursor(input%vars(i))
      end do

      ok = .true.
   end subroutine load_input_data

   subroutine build_full_window(spec, scan)
      type(DataSpec),   intent(inout) :: spec
      type(NetcdfScan), intent(in)    :: scan

      if (allocated(spec%idx_window)) deallocate(spec%idx_window)
      allocate(spec%idx_window(2, 1))

      spec%idx_window(1, 1) = scan%i0
      spec%idx_window(2, 1) = scan%i1
   end subroutine build_full_window


   subroutine ensure_input_data_size(input, nvars)
      type(InputData), intent(inout) :: input
      integer,         intent(in)    :: nvars

      if (.not. allocated(input%vars)) then
         allocate(input%vars(nvars))
      else if (size(input%vars) /= nvars) then
         deallocate(input%vars)
         allocate(input%vars(nvars))
      end if
   end subroutine ensure_input_data_size


   subroutine set_constant_series(series, name, value)
      type(DataVarSeries), intent(inout) :: series
      character(*),        intent(in)    :: name
      real(rk),            intent(in)    :: value

      series%name        = trim(name)
      series%is_const    = .true.
      series%const_value = value
      series%idx         = 1
      series%n           = 0
      series%t_next      = INF_EDGE

      series%time_mode   = DATA_TIME_ABSOLUTE
      series%repeat_year = -1
      series%sim_offset  = 0_lk

      if (allocated(series%t_axis)) deallocate(series%t_axis)
      if (allocated(series%t_edge)) deallocate(series%t_edge)
      if (allocated(series%values)) deallocate(series%values)
   end subroutine set_constant_series


   subroutine init_series_cursor(series)
      type(DataVarSeries), intent(inout) :: series

      if (series%is_const) then
         series%n = 0
         series%idx = 1
         series%t_next = INF_EDGE
         return
      end if

      if (allocated(series%t_axis)) then
         series%n = size(series%t_axis)
      else
         series%n = 0
      end if

      series%idx = 1

      if (series%n <= 1) then
         series%t_next = INF_EDGE
      else
         series%t_next = series%t_edge(2)
      end if
   end subroutine init_series_cursor


   real(rk) function value_at_step(series, model_time) result(value)
      type(DataVarSeries), intent(inout) :: series
      integer(lk),         intent(in)    :: model_time

      if (series%is_const) then
         value = series%const_value
         return
      end if

      if (series%n <= 0) then
         error stop 'value_at_step: empty input series for '//trim(series%name)
      end if

      if (model_time < series%t_edge(series%idx) .or. &
         model_time >= series%t_edge(series%idx + 1)) then
         series%idx = 1
         series%t_next = series%t_edge(2)
      end if

      do while (series%idx < series%n .and. model_time >= series%t_edge(series%idx + 1))
         series%idx = series%idx + 1

         if (series%idx < series%n) then
            series%t_next = series%t_edge(series%idx + 1)
         else
            series%t_next = INF_EDGE
         end if
      end do

      value = series%values(series%idx)
   end function value_at_step


   subroutine print_data_summary(state)
      use, intrinsic :: iso_fortran_env, only: output_unit
      type(DataLoaderState), intent(in) :: state

      integer :: i
      character(:), allocatable :: list_file, list_const, list_off, list_compute
      logical :: any_file, any_const, any_off, any_compute

      list_file    = ''
      list_const   = ''
      list_off     = ''
      list_compute = ''

      any_file    = .false.
      any_const   = .false.
      any_off     = .false.
      any_compute = .false.

      do i = 1, size(state%specs)
         select case (state%specs(i)%input_type)
         case (DATA_INPUT_FILE)
            call append_name(list_file, state%specs(i)%name, any_file)
         case (DATA_INPUT_CONSTANT)
            call append_name(list_const, state%specs(i)%name, any_const)
         case (DATA_INPUT_OFF)
            call append_name(list_off, state%specs(i)%name, any_off)
         case (DATA_INPUT_COMPUTE)
            call append_name(list_compute, state%specs(i)%name, any_compute)
         end select
      end do

      write(output_unit,'(A)') '  ✓ Input data configuration scanned - all checks passed'

      if (any_file)    write(output_unit,'(2X,"Variables from file: ",A)') trim(list_file)
      if (any_const)   write(output_unit,'(2X,"Constant variables: ",A)')  trim(list_const)
      if (any_compute) write(output_unit,'(2X,"Computed variables: ",A)')  trim(list_compute)
      if (any_off)     write(output_unit,'(2X,"Off variables: ",A)')       trim(list_off)

      write(output_unit,'(2X,"Calendar: ",A)') trim(state%sim_cal%name())

      if (state%cfg%load_yearly) then
         write(output_unit,'(2X,"Loading mode: yearly")')
      else
         write(output_unit,'(2X,"Loading mode: full")')
      end if
   end subroutine print_data_summary


   subroutine append_name(list, name, any)
      character(:), allocatable, intent(inout) :: list
      character(*),              intent(in)    :: name
      logical,                   intent(inout) :: any

      if (any) then
         list = trim(list)//', '//trim(name)
      else
         list = trim(name)
         any = .true.
      end if
   end subroutine append_name


   subroutine build_file_list(specs, files)
      type(DataSpec),     intent(in)  :: specs(:)
      type(DataFileInfo), allocatable, intent(out) :: files(:)

      character(:), allocatable :: paths(:)
      integer :: i, nfiles, maxlen
      logical :: found

      maxlen = 1

      do i = 1, size(specs)
         if (specs(i)%input_type == DATA_INPUT_FILE) then
            maxlen = max(maxlen, len_trim(specs(i)%path))
         end if
      end do

      allocate(character(len=maxlen) :: paths(size(specs)))
      paths = ''
      nfiles = 0

      do i = 1, size(specs)
         if (specs(i)%input_type /= DATA_INPUT_FILE) cycle

         found = string_in_list(paths, nfiles, specs(i)%path)

         if (.not. found) then
            nfiles = nfiles + 1
            paths(nfiles) = trim(specs(i)%path)
         end if
      end do

      allocate(files(nfiles))

      do i = 1, nfiles
         files(i)%path      = trim(paths(i))
         files(i)%extension = get_file_extension(paths(i))
         files(i)%format    = infer_format_from_extension(files(i)%extension)
      end do
   end subroutine build_file_list


   logical function string_in_list(list, n, value) result(found)
      character(*), intent(in) :: list(:)
      integer,      intent(in) :: n
      character(*), intent(in) :: value

      integer :: i

      found = .false.

      do i = 1, n
         if (trim(list(i)) == trim(value)) then
            found = .true.
            return
         end if
      end do
   end function string_in_list


   integer function find_file_index(files, path) result(idx)
      type(DataFileInfo), intent(in) :: files(:)
      character(*),       intent(in) :: path

      integer :: i

      idx = 0

      do i = 1, size(files)
         if (trim(files(i)%path) == trim(path)) then
            idx = i
            return
         end if
      end do
   end function find_file_index


   pure integer function find_data_index(input, name) result(idx)
      type(InputData), intent(in) :: input
      character(*),    intent(in) :: name

      integer :: i

      idx = 0
      if (.not. allocated(input%vars)) return

      do i = 1, size(input%vars)
         if (trim(input%vars(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do
   end function find_data_index

   subroutine normalise_specs(specs)
      type(DataSpec), intent(inout) :: specs(:)

      integer :: i

      do i = 1, size(specs)
         if (.not. allocated(specs(i)%name)) specs(i)%name = ''

         select case (specs(i)%input_type)

         case (DATA_INPUT_FILE)
            if (.not. allocated(specs(i)%source_var)) specs(i)%source_var = specs(i)%name
            if (.not. allocated(specs(i)%time_var))   specs(i)%time_var   = 'time'
            if (.not. allocated(specs(i)%path))       specs(i)%path       = ''

            if (specs(i)%format == DATA_FORMAT_UNKNOWN) then
               specs(i)%format = infer_format_from_extension(get_file_extension(specs(i)%path))
            end if

         case (DATA_INPUT_CONSTANT)
            if (.not. allocated(specs(i)%source_var)) specs(i)%source_var = ''
            if (.not. allocated(specs(i)%time_var))   specs(i)%time_var   = ''

         case (DATA_INPUT_COMPUTE, DATA_INPUT_OFF)
            if (.not. allocated(specs(i)%source_var)) specs(i)%source_var = ''
            if (.not. allocated(specs(i)%time_var))   specs(i)%time_var   = ''
            if (.not. allocated(specs(i)%path))       specs(i)%path       = ''

         end select

         if (.not. allocated(specs(i)%units)) specs(i)%units = ''
      end do
   end subroutine normalise_specs


   subroutine validate_specs(specs, ok, errmsg)
      type(DataSpec), intent(in)  :: specs(:)
      logical,        intent(out) :: ok
      character(*),   intent(out) :: errmsg

      integer :: i

      ok = .false.
      errmsg = ''

      do i = 1, size(specs)
         if (.not. allocated(specs(i)%name)) then
            errmsg = 'DataSpec has an unallocated name.'
            return
         end if

         if (len_trim(specs(i)%name) == 0) then
            errmsg = 'DataSpec has an empty name.'
            return
         end if

         select case (specs(i)%input_type)

         case (DATA_INPUT_FILE)
            if (.not. allocated(specs(i)%path) .or. len_trim(specs(i)%path) == 0) then
               errmsg = 'DataSpec '//trim(specs(i)%name)//' has input_type=file but no path.'
               return
            end if

            if (.not. allocated(specs(i)%source_var) .or. len_trim(specs(i)%source_var) == 0) then
               errmsg = 'DataSpec '//trim(specs(i)%name)//' has input_type=file but no source_var.'
               return
            end if

            if (specs(i)%format == DATA_FORMAT_UNKNOWN) then
               errmsg = 'DataSpec '//trim(specs(i)%name)//' has input_type=file but unknown format.'
               return
            end if

         case (DATA_INPUT_CONSTANT)
            ! const_value is already stored.

         case (DATA_INPUT_COMPUTE)
            ! The manager, not the generic loader, should resolve this.

         case (DATA_INPUT_OFF)
            ! Valid. It becomes a constant zero series.

         case default
            errmsg = 'DataSpec '//trim(specs(i)%name)//' has an invalid input_type.'
            return
         end select
      end do

      ok = .true.
   end subroutine validate_specs


   function get_file_extension(path) result(ext)
      character(*), intent(in) :: path
      character(:), allocatable :: ext

      integer :: i, lastdot

      lastdot = 0

      do i = len_trim(path), 1, -1
         if (path(i:i) == '.') then
            lastdot = i
            exit
         end if
      end do

      if (lastdot == 0 .or. lastdot == len_trim(path)) then
         ext = ''
      else
         ext = to_lower(path(lastdot+1:len_trim(path)))
      end if

      if (ext == 'nc4' .or. ext == 'cdf') ext = 'nc'
   end function get_file_extension


   integer function infer_format_from_extension(extension) result(format)
      character(*), intent(in) :: extension

      select case (trim(to_lower(extension)))
      case ('nc', 'netcdf')
         format = DATA_FORMAT_NETCDF
      case ('csv')
         format = DATA_FORMAT_CSV
      case default
         format = DATA_FORMAT_UNKNOWN
      end select
   end function infer_format_from_extension


   subroutine clear_state(state)
      type(DataLoaderState), intent(inout) :: state

      state%cfg%repeat_enabled = .false.
      state%cfg%repeat_year    = -huge(1)
      state%cfg%cfg_calendar   = cal_unknown
      state%cfg%load_yearly    = .true.

      state%sim_y_start = 0
      state%sim_y_end   = 0
      state%sim_cal%kind = cal_unknown

      if (allocated(state%specs))    deallocate(state%specs)
      if (allocated(state%files))    deallocate(state%files)
      if (allocated(state%nc_scans)) deallocate(state%nc_scans)
   end subroutine clear_state

end module data_loader