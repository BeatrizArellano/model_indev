! Useful functions for time
module time_utils
  use, intrinsic :: iso_fortran_env, only: error_unit
  use precision_types, only: rk, ik                  ! Importing real64 and int32
  use time_types,      only: DateTime, cal_gregorian, cal_proleptic, cal_noleap, cal_all_leap, cal_360_day
  use str_utils,       only: replace_char_inplace
  use stats_utils,     only: median_rk

  implicit none
  private




  public :: is_leap_gregorian, check_time_monotonic
  public :: parse_datetime_str, datetime_from_string, datetime_to_str
  public :: validate_datetime, validate_calendar
  public :: is_datetime_before, is_datetime_equal, is_datetime_after
  public :: index_at_or_before, index_at_or_after, index_linear_weights, detect_frequency

contains

    ! Checks if a given year is leap in the Gregorian calendar
    pure logical function is_leap_gregorian(year) result(yes)
        integer, intent(in) :: year
        yes = (mod(year,4)==0 .and. (mod(year,100)/=0 .or. mod(year,400)==0))
    end function is_leap_gregorian

    ! Checks if datetime a is before datetime b
    pure elemental logical function is_datetime_before(a, b) result(ans)
        type(DateTime), intent(in) :: a, b
        ans = .false.
        if (a%year  /= b%year ) then; ans = a%year  < b%year;  return; end if
        if (a%month /= b%month) then; ans = a%month < b%month; return; end if
        if (a%day   /= b%day  ) then; ans = a%day   < b%day;   return; end if
        if (a%hour  /= b%hour ) then; ans = a%hour  < b%hour;  return; end if
        if (a%minute/= b%minute)then; ans = a%minute< b%minute;return; end if
        if (a%second/= b%second)then; ans = a%second< b%second;return; end if
        ! equal down to seconds -> not before
    end function is_datetime_before

    ! Checks if two date-times are equal. 
    ! If include_has_time=.true., equality also requires a%has_time==b%has_time.
    ! Default (absent or .false.) ignores has_time (i.e., equality by instant).
    pure elemental logical function is_datetime_equal(a, b, include_has_time) result(ans)
        type(DateTime), intent(in) :: a, b
        logical,       intent(in), optional :: include_has_time
        logical :: same_instant, strict

        same_instant =  (a%year   == b%year  ) .and. &
                        (a%month  == b%month ) .and. &
                        (a%day    == b%day   ) .and. &
                        (a%hour   == b%hour  ) .and. &
                        (a%minute == b%minute) .and. &
                        (a%second == b%second)

        strict = present(include_has_time) .and. include_has_time
        if (strict) then
            ans = same_instant .and. (a%has_time .eqv. b%has_time)
        else
            ans = same_instant
        end if
    end function is_datetime_equal


    pure elemental logical function is_datetime_after(a, b) result(ans)
        type(DateTime), intent(in) :: a, b
        ans = is_datetime_before(b, a)
    end function is_datetime_after


    ! Parses ISO-ish date-time string and returns year, month, day, hour
    ! Accepts "YYYY-MM-DD", "YYYY-MM-DD hh:mm", "YYYY-MM-DD hh:mm:ss",
    ! "YYYY-MM-DDThh:mm[:ss]" and optional trailing 'Z' (UTC) for ISO 8601 dates
    ! On success: returns Y,M,D,H,Mn,S with missing time parts = 0; ok=.true.
    ! If time is present then has_time=true
    ! On failure: ok=false and errmsg describing the problem.
    subroutine parse_datetime_str(s_in, year, month, day, hour, minute, second, has_time, ok, errmsg)
        character(*), intent(in)  :: s_in
        integer,      intent(out) :: year, month, day, hour, minute, second
        logical,      intent(out) :: has_time, ok
        character(*), intent(out) :: errmsg

        character(len=:), allocatable :: s, datepart, timepart, buf
        integer :: lt, sep, ios

        year=0; month=0; day=0; hour=0; minute=0; second=0
        ok=.false.; has_time=.false.; errmsg=''

        ! Normalize: trim, strip quotes, make 'T' a space, strip trailing 'Z'
        s = trim(adjustl(s_in))
        lt = len_trim(s)
        if (lt>=2) then
            if ( (s(1:1)=='"'  .and. s(lt:lt)=='"')  .or. &
                (s(1:1)=='''' .and. s(lt:lt)=='''') ) then
            s = s(2:lt-1)
            end if
        end if
        call replace_char_inplace(s, 'T', ' ')
        lt = len_trim(s)
        if (lt>0 .and. s(lt:lt)=='Z') s = s(:lt-1)
        s = trim(s)

        ! Find first whitespace (space or tab) splitting date and time
        sep = 0
        do lt=1,len(s)
            if (s(lt:lt)==' ' .or. s(lt:lt)==char(9)) then
            sep = lt
            exit
            end if
        end do

        if (sep>0) then
            datepart = trim(s(:sep-1))
            timepart = trim(adjustl(s(sep+1:)))
        else
            datepart = s
            timepart = ''
        end if

        ! Sanitize: map any non-digit to a single space (robust to weird dashes/NBSP)
        call non_digits_to_space(datepart, buf); datepart = trim(buf)
        if (len(timepart)>0) then
            call non_digits_to_space(timepart, buf); timepart = trim(buf)
        end if

        ! Read date as three integers with list-directed input
        read(datepart,*,iostat=ios) year, month, day
        if (ios /= 0) then
            errmsg = 'CF time: could not parse date part "'//trim(s_in)//'".'
            return
        end if

        if (len_trim(timepart) > 0) then
            ! Try HH MM SS, then HH MM
            read(timepart,*,iostat=ios) hour, minute, second
            if (ios /= 0) then
            second = 0
            read(timepart,*,iostat=ios) hour, minute
            if (ios /= 0) then
                errmsg = 'CF time: could not parse time part "'//trim(s_in)//'".'
                return
            end if
            end if
            has_time = .true.; ok = .true.; return
        else
            has_time = .false.; ok = .true.; return
        end if

        contains
        pure subroutine non_digits_to_space(a, out)
            character(*), intent(in)  :: a
            character(len=:), allocatable, intent(out) :: out
            integer :: i, n
            logical :: last_space
            character(len=:), allocatable :: tmp

            n = len_trim(a)
            if (n==0) then
            out = ''; return
            end if
            tmp = repeat(' ', n*2)  ! generous buffer
            last_space = .true.; i = 0

            block
            integer :: k, ich
            do k=1,n
                ich = iachar(a(k:k))
                if (ich>=iachar('0') .and. ich<=iachar('9')) then
                i = i + 1
                tmp(i:i) = a(k:k)
                last_space = .false.
                else
                if (.not. last_space) then
                    i = i + 1
                    tmp(i:i) = ' '
                    last_space = .true.
                end if
                end if
            end do
            end block

            if (i<=0) then
            out = ''
            else
            out = trim(tmp(:i))
            end if
        end subroutine non_digits_to_space
    end subroutine parse_datetime_str
    

    ! Validates whether a date-time is correct after parsing the date-time string
    pure subroutine validate_datetime(y,m,d,h,mi,s, ok, errmsg, calendar)
        integer, intent(in) :: y,m,d,h,mi,s
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg
        integer, intent(in), optional :: calendar
        integer :: cal, daysinmonth

        cal = merge(calendar, cal_gregorian, present(calendar))  ! default

        ok = .false.; errmsg=''
        
        if (m < 1 .or. m > 12)           then; errmsg='Invalid month [1..12]';       return; end if

        daysinmonth = days_in_month_calendar(y,m,cal)
        if (daysinmonth <= 0)            then; errmsg='Invalid number of days';      return; end if
        if (d < 1 .or. d > daysinmonth)  then; errmsg='Invalid day for month using that calendar.'; return; end if
        if (h < 0 .or. h > 23)           then; errmsg='Invalid hour [0..23]';        return; end if
        if (mi < 0 .or. mi > 59)         then; errmsg='Invalid minute [0..59]';      return; end if
        if (s < 0 .or. s > 59)           then; errmsg='Invalid second [0..59]';      return; end if

        ok = .true.
    end subroutine validate_datetime

    pure subroutine validate_calendar(cal_code, ok)
        integer, intent(in)  :: cal_code
        logical, intent(out) :: ok
        select case (cal_code)
        case (cal_gregorian, cal_proleptic, cal_noleap, cal_all_leap, cal_360_day)
            ok = .true.
        case default
            ok = .false.
        end select
    end subroutine validate_calendar


    pure integer function days_in_month_calendar(y,m,cal) result(d)
        integer, intent(in) :: y,m,cal
        integer, parameter :: dm_greg(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
        select case (cal)
        case (cal_gregorian,cal_proleptic)
            d = dm_greg(m); if (m==2 .and. is_leap_gregorian(y)) d = 29
        case (cal_noleap)
            d = dm_greg(m)                                      ! Feb=28 always
        case (cal_all_leap)
            d = dm_greg(m); if (m==2) d = 29                    ! Feb=29 always
        case (cal_360_day)
            d = 30                                              ! 12x30 months
        case default
            d = -1                                              ! unknown calendar
        end select
    end function days_in_month_calendar

    ! Wrapper for parse_datetime_str to stop if there is an error
    ! Otherwise assign the parsed values to a DateTime structure
    function datetime_from_string(datetime_str, ok, errmsg, calendar) result(dati)
        character(*), intent(in)  :: datetime_str
        logical,      intent(out) :: ok
        character(*), intent(out) :: errmsg
        integer,      intent(in),  optional :: calendar
        type(DateTime) :: dati

        integer :: y, mon, d, h, mi, sec
        logical :: has_time, okp
        character(len=256) :: msg
        character(len=:), allocatable :: stopmsg

        ! Parsing
        call parse_datetime_str(datetime_str, y, mon, d, h, mi, sec, has_time, okp, msg)
        if (.not. okp) then
            ok = .false.; errmsg = msg
            stopmsg = 'Error parsing datetime: "'//trim(datetime_str)//'" - '//trim(msg)
            error stop stopmsg
        end if

        ! Validation - Stop if the datetime is invalid
        call validate_datetime(y, mon, d, h, mi, sec, okp, msg, calendar)
        if (.not. okp) then
            ok = .false.; errmsg = msg
            stopmsg = 'Datetime "'//trim(datetime_str)//'" is not valid: '//trim(msg)
            error stop stopmsg
        end if

        ! Assign data
        dati%year   = y
        dati%month  = mon
        dati%day    = d
        dati%hour   = h
        dati%minute = mi
        dati%second = sec
        dati%has_time = has_time
        ok = .true.; errmsg = ''
    end function datetime_from_string

    !> Convert a DateTime to "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" (if has_time=.true.)
    pure function datetime_to_str(dt) result(s)
        type(DateTime), intent(in) :: dt
        character(:), allocatable  :: s
        character(len=19) :: buf19
        character(len=10) :: buf10

        if (dt%has_time) then
        write(buf19,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
                dt%year, dt%month, dt%day, dt%hour, dt%minute, dt%second
        s = trim(buf19)
        else
        write(buf10,'(I4.4,"-",I2.2,"-",I2.2)') dt%year, dt%month, dt%day
        s = trim(buf10)
        end if
    end function datetime_to_str


    !-----------------------------------------------------------
    ! Returns index in time array of requested time (in seconds) 
    !   t_s: time-array in seconds
    !   t_query: time in seconds to find in array
    ! Returns i such that t_s(i) <= t_query < t_s(i+1) -> Earliest index found
    !     IMPORTANT: If t_query ≤ t_s(1) returns 1
    !                If t_query ≥ t_s(n) returns n
    !-----------------------------------------------------------
    integer function index_at_or_before(t_s, t_query) result(i0)
        real(rk), intent(in) :: t_s(:), t_query
        integer :: lo, hi, mid, n
        n = size(t_s)
        if (n == 0) then
            i0 = 0; return
        end if
        if (t_query <= t_s(1)) then
            i0 = 1; return
        end if
        if (t_query >= t_s(n)) then
            i0 = n; return
        end if
        lo = 1; hi = n
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (t_s(mid) <= t_query) then
                lo = mid
            else
                hi = mid
            end if
        end do
        i0 = lo
    end function index_at_or_before

    pure integer function index_at_or_after(t_s, t_query) result(i1)
        real(rk), intent(in) :: t_s(:)
        real(rk), intent(in) :: t_query
        integer :: lo, hi, mid, n

        n = size(t_s)
        if (n == 0) then
            i1 = 0             ! empty input
            return
        end if
        if (t_query <= t_s(1)) then
            i1 = 1
            return
        end if
        if (t_query > t_s(n)) then
            i1 = n + 1         ! sentinel: "one past the last"
            return
        end if

        ! Binary search for first i with t_s(i) >= t_query
        lo = 1; hi = n
        do while (lo < hi)
            mid = (lo + hi) / 2
            if (t_s(mid) < t_query) then
            lo = mid + 1
            else
            hi = mid
            end if
        end do
        i1 = lo
    end function index_at_or_after


    !----------------------------------------------------------------------------
    ! Given a time (t_query), index_linear_weights finds bracketing indices 
    ! and linear interpolation weights for t_query on an ascending time array (t_s).
    !   t_s: time-array in seconds, t_query: time in seconds to find in array
    ! Returns:
    ! i0, i1: Bracketing indices
    ! w0, w1: weights for those indices to be used as:
    !       y ≈ w0 * y(i0) + w1 * y(i1)
    !-----------------------------------------------------------------------
    subroutine index_linear_weights(t_s, t_query, i0, i1, w0, w1)
        real(rk), intent(in)  :: t_s(:), t_query
        integer,  intent(out) :: i0, i1
        real(rk), intent(out) :: w0, w1
        integer :: n
        real(rk) :: t0, t1, denom

        n = size(t_s)
        if (n == 0) then
            i0 = 0; i1 = 0; w0 = 0.0_rk; w1 = 0.0_rk; return
        end if
        if (t_query <= t_s(1)) then
            i0 = 1; i1 = 1; w0 = 1.0_rk; w1 = 0.0_rk; return
        end if
        if (t_query >= t_s(n)) then
            i0 = n; i1 = n; w0 = 1.0_rk; w1 = 0.0_rk; return
        end if

        i0 = index_at_or_before(t_s, t_query)
        i1 = i0 + 1
        t0 = t_s(i0); t1 = t_s(i1)
        denom = max(t1 - t0, 1.0e-12_rk)
        w1 = (t_query - t0) / denom
        w0 = 1.0_rk - w1
    end subroutine index_linear_weights

    ! Checks whether a time-series is monotonic
    subroutine check_time_monotonic(t_s, is_monotonic, has_equal_consecutive, i_bad)
        use precision_types, only: rk
        real(rk), intent(in)  :: t_s(:)
        logical,  intent(out) :: is_monotonic, has_equal_consecutive
        integer,  intent(out) :: i_bad

        integer :: i
        real(rk) :: d
        real(rk), parameter :: tol_neg = -1.0e-6_rk  ! allow tiny noise; anything < tol_neg is an inversion
        real(rk), parameter :: tol_tie =  1.0e-9_rk  ! treat <= 1 ns as a "tie"

        is_monotonic = .true.
        has_equal_consecutive     = .false.
        i_bad        = 0

        do i = 1, size(t_s)-1
            d = t_s(i+1) - t_s(i)
            if (d < tol_neg) then
                is_monotonic = .false.
                i_bad = i
                return
            end if
            if (abs(d) <= tol_tie) has_equal_consecutive = .true.
        end do
    end subroutine check_time_monotonic


    !---------------------------------------------------------------------------------------------
    ! Infers the time step of a time array and checks regularity: 
    !             step_s = median(Δt), is_regular = all(|Δt−step_s| ≤ tolerance 0.01 s).
    ! Returns the Relative maximum deviation too:
    ! rel_max_dev = max(|Δt−step_s|)/step_s 
    !---------------------------------------------------------------------------------------------
    subroutine detect_frequency(t_s, step_s, is_regular, rel_max_dev)
        real(rk), intent(in)  :: t_s(:)     ! seconds since reference date
        real(rk), intent(out) :: step_s     ! canonical step in seconds (median diff)
        logical,  intent(out) :: is_regular
        real(rk), intent(out) :: rel_max_dev ! Relative maximum deviation (max |diff - step| / step)

        integer :: n, i
        real(rk), allocatable :: d(:)
        real(rk) :: med, max_dev, tol_abs

        n = size(t_s)
        if (n < 2) then
            step_s = 0.0_rk; is_regular = .true.; rel_max_dev = 0.0_rk
            return
        end if
        allocate(d(n-1))
        do i=1, n-1
            d(i) = t_s(i+1) - t_s(i)
        end do

        med = median_rk(d)
        step_s = med

        max_dev = 0.0_rk
        do i=1, size(d)
            max_dev = max(max_dev, abs(d(i) - step_s))
        end do
        tol_abs   = max(1.0e-6_rk*max(1.0_rk,step_s), 0.01_rk)  ! 1e-6 relative or 0.01 s
        is_regular = all(abs(d - step_s) <= tol_abs)
        if (step_s > 0.0_rk) then
            rel_max_dev = max_dev / step_s
        else
            rel_max_dev = 0.0_rk
        end if
    end subroutine detect_frequency


end module time_utils