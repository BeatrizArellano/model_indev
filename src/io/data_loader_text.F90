module data_loader_text
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   use cf_time_utils,   only: seconds_between_datetimes
   use data_types,      only: DataSpec, DataVarSeries, DATA_INPUT_FILE
   use precision_types, only: rk, lk
   use text_parser,     only: TextTable, parse_text_table
   use time_types,      only: DateTime, CFCalendar, CFUnits, TimeAxis
   use time_utils,      only: parse_datetime_str, validate_datetime, &
                              check_time_monotonic, detect_frequency, &
                              index_at_or_before

   implicit none
   private

   public :: TextScan
   public :: scan_text_file
   public :: load_text_series
   public :: build_text_full_window

   integer(lk), parameter :: INF_EDGE = huge(1_lk)

   type :: TextScan
      character(:), allocatable :: path
      character(:), allocatable :: time_name

      type(CFCalendar) :: cal
      type(CFUnits)    :: u
      type(TimeAxis)   :: axis

      type(DateTime), allocatable :: datetimes(:)

      integer :: i0 = 1
      integer :: i1 = 1

      real(rk) :: median_dt = 0.0_rk
      logical  :: is_regular = .true.
      real(rk) :: rel_max_dev = 0.0_rk
   end type TextScan

contains

   subroutine scan_text_file(path, specs, calendar_cfg, start_datetime, end_datetime, scan, ok, errmsg)
        character(*),    intent(in)  :: path
        type(DataSpec),  intent(in)  :: specs(:)
        type(CFCalendar), intent(in) :: calendar_cfg
        type(DateTime),  intent(in)  :: start_datetime, end_datetime
        type(TextScan),  intent(out) :: scan
        logical,         intent(out) :: ok
        character(*),    intent(out) :: errmsg

        type(TextTable) :: table
        integer :: itime, i, ntime, i_bad
        logical :: lok, has_equal
        real(rk) :: t_start, t_end

        ok = .false.
        errmsg = ''
        call clear_text_scan(scan)

        scan%path = trim(path)
        scan%cal  = calendar_cfg
        call set_text_units_to_sim_start(scan%u, start_datetime)
        scan%axis%cal = scan%cal
        scan%axis%u   = scan%u

        call parse_text_table(path, table)

        if (table%find_column('lat') > 0 .or. table%find_column('latitude') > 0 .or. &
            table%find_column('lon') > 0 .or. table%find_column('longitude') > 0) then
            write(*,'(A)') 'WARNING DataLoader: Text file '//trim(path)// &
                  ' contains latitude/longitude columns, but plain text input is treated as a single-point time series; lat/lon columns are ignored.'
        end if

        call find_text_time_column(table, specs, path, scan%time_name, itime, ok, errmsg)
        if (.not. ok) then
            call table%free()
            return
        end if

        ntime = table%nrow
        if (ntime <= 0) then
            errmsg = 'Text file '//trim(path)//' contains no data rows.'
            call table%free()
            return
        end if

        allocate(scan%datetimes(ntime))
        allocate(scan%axis%t_s(ntime))

        do i = 1, ntime
            call parse_text_datetime(trim(table%values(i, itime)), scan%cal, scan%datetimes(i), lok, errmsg)
            if (.not. lok) then
                errmsg = 'Invalid time in '//trim(path)//' at data row '//trim(itoa(i))//': '//trim(errmsg)
                call table%free()
                return
            end if

            scan%axis%t_s(i) = real(seconds_between_datetimes(scan%cal, scan%datetimes(i), start_datetime), rk)
        end do

        if (any(.not. ieee_is_finite(scan%axis%t_s))) then
            errmsg = 'Text file '//trim(path)//' contains invalid time values.'
            call table%free()
            return
        end if

        call check_time_monotonic(scan%axis%t_s, lok, has_equal, i_bad)
        if (.not. lok) then
            errmsg = 'Text time is not monotonic in '//trim(path)//' near rows '// &
                    trim(itoa(i_bad))//' and '//trim(itoa(i_bad + 1))//'.'
            call table%free()
            return
        end if

        if (has_equal) then
            errmsg = 'Text file '//trim(path)//' contains duplicate consecutive timestamps.'
            call table%free()
            return
        end if

        call detect_frequency(scan%axis%t_s, scan%median_dt, scan%is_regular, scan%rel_max_dev)

        scan%axis%t_first = scan%axis%t_s(1)
        scan%axis%t_last  = scan%axis%t_s(ntime)

        t_start = 0.0_rk
        t_end   = real(seconds_between_datetimes(scan%cal, end_datetime, start_datetime), rk)

        if (scan%median_dt <= 0.0_rk) then
            errmsg = 'Detected non-positive text time step in '//trim(path)//'.'
            call table%free()
            return
        end if

        if (t_start < scan%axis%t_first) then
            errmsg = 'Text file '//trim(path)//' starts after the requested simulation period.'
            call table%free()
            return
        end if

        if (t_end > scan%axis%t_last + scan%median_dt) then
            errmsg = 'Text file '//trim(path)//' ends before the requested simulation period is covered.'
            call table%free()
            return
        end if

        scan%i0 = max(1, index_at_or_before(scan%axis%t_s, t_start))
        scan%i1 = max(scan%i0, index_at_or_before(scan%axis%t_s, t_end))

        call validate_text_columns_and_values(table, specs, path, scan%i0, scan%i1, ok, errmsg)
        if (.not. ok) then
            call table%free()
            return
        end if

        call table%free()
        ok = .true.
    end subroutine scan_text_file


    subroutine build_text_full_window(spec, scan)
        type(DataSpec), intent(inout) :: spec
        type(TextScan), intent(in)    :: scan

        if (allocated(spec%idx_window)) deallocate(spec%idx_window)
        allocate(spec%idx_window(2, 1))

        spec%idx_window(1, 1) = scan%i0
        spec%idx_window(2, 1) = scan%i1
    end subroutine build_text_full_window


    subroutine load_text_series(series, scan, spec, i0, i1)
        type(DataVarSeries), intent(inout) :: series
        type(TextScan),      intent(in)    :: scan
        type(DataSpec),      intent(in)    :: spec
        integer,             intent(in)    :: i0, i1

        type(TextTable) :: table
        integer     :: jcol, i, nt, ios
        integer(lk) :: dt_last
        character(len=512) :: msg

        nt = max(0, i1 - i0 + 1)
        if (nt <= 0) error stop 'load_text_series: empty time window for '//trim(spec%name)

        call parse_text_table(scan%path, table)

        jcol = table%find_column(trim(spec%source_var))
        if (jcol <= 0) then
            call table%free()
            error stop 'load_text_series: missing column '//trim(spec%source_var)
        end if

        call clear_series_arrays(series)

        allocate(series%t_axis(nt))
        allocate(series%t_edge(nt + 1))
        allocate(series%values(nt))

        series%name     = trim(spec%name)
        series%units    = trim(spec%units)
        series%is_const = .false.
        series%n        = nt
        series%idx      = 1

        series%cal         = scan%cal
        series%u           = scan%u
        series%sim_offset  = 0_lk
        series%time_mode   = spec%time_mode
        series%repeat_year = spec%repeat_year

        series%t_axis = int(nint(scan%axis%t_s(i0:i1)), lk)

        do i = 1, nt
            read(table%values(i0 + i - 1, jcol), *, iostat=ios) series%values(i)
            if (ios /= 0 .or. .not. ieee_is_finite(series%values(i))) then
                write(msg,'(A,A,A,I0)') 'load_text_series: invalid numeric value in column ', &
                                        trim(spec%source_var), ' at data row ', i0 + i - 1
                call table%free()
                error stop trim(msg)
            end if
        end do

        series%t_edge(1:nt) = series%t_axis(:)

        if (nt >= 2) then
            dt_last = max(1_lk, series%t_axis(nt) - series%t_axis(nt - 1))
            series%t_next = series%t_edge(2)
        else
            dt_last = max(1_lk, int(nint(scan%median_dt), lk))
            series%t_next = INF_EDGE
        end if

        series%t_edge(nt + 1) = series%t_axis(nt) + dt_last

        call table%free()
        return
    end subroutine load_text_series


    subroutine parse_text_datetime(txt, cal, dt, ok, errmsg)
        character(*),     intent(in)  :: txt
        type(CFCalendar), intent(in)  :: cal
        type(DateTime),   intent(out) :: dt
        logical,          intent(out) :: ok
        character(*),     intent(out) :: errmsg

        call parse_datetime_str(trim(txt), dt%year, dt%month, dt%day, &
                                dt%hour, dt%minute, dt%second, &
                                dt%has_time, ok, errmsg)
        if (.not. ok) return

        call validate_datetime(dt%year, dt%month, dt%day, &
                                dt%hour, dt%minute, dt%second, &
                                ok, errmsg, cal%kind)
    end subroutine parse_text_datetime


    subroutine set_text_units_to_sim_start(u, start_datetime)
        type(CFUnits),   intent(out) :: u
        type(DateTime),  intent(in)  :: start_datetime

        u%timeunit_to_seconds = 1.0_rk
        u%reference_year      = start_datetime%year
        u%reference_month     = start_datetime%month
        u%reference_day       = start_datetime%day
        u%reference_hour      = start_datetime%hour
        u%reference_min       = start_datetime%minute
        u%reference_sec       = start_datetime%second
        u%has_time            = .true.
    end subroutine set_text_units_to_sim_start

   !================================================================
   !   Internal
   !================================================================

    subroutine find_text_time_column(table, specs, path, time_name, itime, ok, errmsg)
        type(TextTable), intent(in) :: table
        type(DataSpec),  intent(in) :: specs(:)
        character(*),    intent(in) :: path
        character(:), allocatable, intent(out) :: time_name
        integer, intent(out) :: itime
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        integer :: i

        ok = .false.
        errmsg = ''
        itime = 0
        time_name = 'time'

        do i = 1, size(specs)
            if (specs(i)%input_type /= DATA_INPUT_FILE) cycle
            if (.not. allocated(specs(i)%path)) cycle
            if (trim(specs(i)%path) /= trim(path)) cycle
            if (allocated(specs(i)%time_var)) then
                if (len_trim(specs(i)%time_var) > 0) then
                    time_name = trim(specs(i)%time_var)
                    exit
                end if
            end if
        end do

        itime = table%find_column(trim(time_name))

        if (itime <= 0 .and. trim(time_name) /= 'time') then
            itime = table%find_column('time')
            if (itime > 0) time_name = 'time'
        end if

        if (itime <= 0) then
            errmsg = 'Text file is missing required time column "'//trim(time_name)//'".'
            return
        end if

        ok = .true.
    end subroutine find_text_time_column

    subroutine validate_text_columns_and_values(table, specs, path, i0, i1, ok, errmsg)
        type(TextTable), intent(in) :: table
        type(DataSpec),  intent(in) :: specs(:)
        character(*),    intent(in) :: path
        integer, intent(in) :: i0, i1
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        integer :: i, j, jcol, ios
        real(rk) :: x

        ok = .false.
        errmsg = ''

        do j = 1, size(specs)
            if (specs(j)%input_type /= DATA_INPUT_FILE) cycle
            if (.not. allocated(specs(j)%path)) cycle
            if (trim(specs(j)%path) /= trim(path)) cycle
            if (.not. allocated(specs(j)%source_var)) cycle
            if (len_trim(specs(j)%source_var) == 0) cycle

            jcol = table%find_column(trim(specs(j)%source_var))
            if (jcol <= 0) then
                errmsg = 'Text file is missing required variable column "'//trim(specs(j)%source_var)//'".'
                return
            end if

            do i = i0, i1
                read(table%values(i, jcol), *, iostat=ios) x
                if (ios /= 0 .or. .not. ieee_is_finite(x)) then
                    errmsg = 'Invalid numeric value for "'//trim(specs(j)%source_var)// &
                                '" at data row '//trim(itoa(i))//'.'
                    return
                end if
            end do
        end do

        ok = .true.
    end subroutine validate_text_columns_and_values

    subroutine clear_text_scan(scan)
        type(TextScan), intent(inout) :: scan

        if (allocated(scan%path)) deallocate(scan%path)
        if (allocated(scan%time_name)) deallocate(scan%time_name)
        if (allocated(scan%datetimes)) deallocate(scan%datetimes)
        if (allocated(scan%axis%t_s)) deallocate(scan%axis%t_s)

        scan%i0 = 1
        scan%i1 = 1
        scan%median_dt = 0.0_rk
        scan%is_regular = .true.
        scan%rel_max_dev = 0.0_rk
    end subroutine clear_text_scan

    subroutine clear_series_arrays(series)
        type(DataVarSeries), intent(inout) :: series

        if (allocated(series%t_axis)) deallocate(series%t_axis)
        if (allocated(series%t_edge)) deallocate(series%t_edge)
        if (allocated(series%values)) deallocate(series%values)
    end subroutine clear_series_arrays

    pure function itoa(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(len=32) :: buf

        write(buf, '(I0)') i
        s = trim(buf)
    end function itoa

end module data_loader_text

