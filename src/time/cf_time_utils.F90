! src/utils/cf_time_utils.F90
! Module containing functions and utilities to parse time represented 
! according to the Climate and Forecast (CF) Metadata Conventions.
!
! This module supports only 5 calendars:
!   Gregorian or standard, Proleptic Greforian,
!   No Leap (always 365 days), All leap (always 366 days), 
!   360 day: 12 months of 30 days each. 
!
! This module can read time units following the CF Metadata convention:
!           "<units> since YYYY-MM-DD hh:mm:ss"
! where units can be seconds, minutes, hours, days, weeks
module cf_time_utils
  use precision_types, only: rk, lk     ! Importing real64 and int64
  use str_utils,       only: to_lower
  use time_utils,      only: is_leap_gregorian, parse_datetime_str                             
  use stats_utils,     only: sort_real_inplace
  use calendar_types,  only: CFCalendar, cal_gregorian, cal_proleptic, cal_noleap, cal_all_leap, cal_360_day


  implicit none
  private

  type, public :: CFUnits
    ! Holds the CF units since the reference point in time (date)
     real(rk) :: timeunit_to_seconds = 0.0_rk    !conversion factor from the declared CF time unit to seconds (e.g. 86400 for "days since")
     integer  :: reference_year = 1990, reference_month = 1, reference_day = 1
     integer  :: reference_hour = 0,    reference_min = 0, reference_sec = 0
     logical  :: has_time = .false.
  end type

  type, public :: TimeAxis
    ! Data structure to hold a time coordinate using seconds as reference date
     type(CFCalendar) :: cal                        ! calendar
     type(CFUnits)    :: u                          ! units
     real(rk), allocatable :: t_s(:)                ! seconds since reference date
     real(rk) :: t_first = 0.0_rk, t_last = 0.0_rk  ! earliest and last time values in in seconds since the reference date
  end type

  public :: parse_cf_time, seconds_since_datetime, seconds_to_datetime
  public :: get_days_in_year  

  ! Cumulative days before each month for a non-leap year and for a leap year
  !Month:      Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
  !cum_nonleap: 0   31   59   90  120  151  181  212  243  273  304  334
  !cum_leap:    0   31   60   91  121  152  182  213  244  274  305  335
  integer, parameter :: cum_nonleap(12) = [0,31,59,90,120,151,181,212,243,273,304,334]
  integer, parameter :: cum_leap(12) = [0,31,60,91,121,152,182,213,244,274,305,335]

contains

    !==================================================================
    ! Parse CF-compliant time units strings: "<unit> since <reference date>"
    ! Determines the calendar
    !================================================================
    subroutine parse_cf_time(units_attr, calendar_attr, default_calendar, u, cal, ok, errmsg)
        character(*), intent(in)  :: units_attr
        character(*), intent(in)  :: calendar_attr
        character(*), intent(in)  :: default_calendar   ! e.g. "gregorian" if attribute missing
        type(CFUnits), intent(out) :: u
        type(CFCalendar), intent(out) :: cal
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        character(len=:), allocatable :: time_unit, reference_date_str, cal_str
        logical :: ok_units, ok_reference_date, ok_cal

        ok = .false.; errmsg = ''
        
        ! Split the units string and returns units (e.g. days) and the reference date string
        call split_units_since(units_attr, time_unit, reference_date_str, ok_units)
        if (.not. ok_units) then
            errmsg = 'CF time: invalid "units" attribute. Expected "<unit> since <YYYY-MM-DD[ [T]hh:mm[:ss][Z]]>".'
            return
        end if
        ! Verifies the time-units and assigns number of seconds per time-unit 
        call verify_time_unit(time_unit, u%timeunit_to_seconds, ok_units)
        if (.not. ok_units) then
            errmsg = 'CF time: unsupported unit "'//trim(time_unit)//'". Accepted: seconds, minutes, hours, days, weeks.'
            return
        end if

        ! Parse the reference date string and fills the corresponding units in u (CFUnits type)
        call parse_datetime_str(reference_date_str, u%reference_year, u%reference_month, u%reference_day, u%reference_hour, &
                                u%reference_min, u%reference_sec, u%has_time, ok_reference_date, errmsg)
        if (.not. ok_reference_date) return

        if (len_trim(calendar_attr) == 0 .and. &
           (len_trim(default_calendar) == 0 .or. to_lower(trim(default_calendar)) == 'unknown' ) ) then
            ok     = .false.
            errmsg = 'Forcing data: missing calendar (no file attribute and parameter in main.yaml). ' //  &
                     'Please add CF ''calendar'' to the file or set time.calendar in main config file.'
            return
        end if

        if (len_trim(calendar_attr) > 0) then
            cal_str = to_lower(trim(calendar_attr))
        else         
            if ( to_lower(trim(default_calendar)) == 'unknown' ) then
                cal_str = ''                            ! treats "unknown" as missing so verify_calendar_str won't see it
            else
                cal_str = to_lower(trim(default_calendar))  ! uses default calendar if calendar_attr=''
            end if       
        end if
        ! Verifies support for the calendar
        call verify_calendar_str(cal_str, cal%kind, ok_cal)
        if (.not. ok_cal) then
            errmsg = 'CF time: unsupported calendar "'//trim(cal_str)//'".'
            return
        end if

        ok = .true.
    end subroutine parse_cf_time

    !-----------------------------------------------------
    ! Convert a date/time to seconds since CF reference date
    !-------------------------------------------------
    pure function seconds_since_datetime(cal, u, year, month, day, hour, minute, second) result(tsec)
        type(CFCalendar), intent(in) :: cal
        type(CFUnits),    intent(in) :: u
        integer, intent(in) :: year, month, day, hour, minute, second
        real(rk) :: tsec                   ! Total number of seconds since reference date-time
        integer(lk) :: numdays_in_date, numdays_in_reference_date, numdays_since_ref
        integer(lk) :: nsec_in_target_time, nsec_in_reftime, delta_sec

        numdays_in_date = date_to_number_of_days(cal%kind, year, month, day)          ! Number of days in target date (since 0001-01-01)
        numdays_in_reference_date  = date_to_number_of_days(cal%kind, u%reference_year, u%reference_month, u%reference_day) ! Number of days in reference date
        numdays_since_ref    = numdays_in_date - numdays_in_reference_date                   ! Number of days since seference date

        nsec_in_target_time = int(hour, lk)*3600_lk + int(minute, lk)*60_lk + int(second, lk)
        nsec_in_reftime  = int(u%reference_hour, lk)*3600_lk + int(u%reference_min, lk)*60_lk + int(u%reference_sec, lk)
        delta_sec    = numdays_since_ref*86400_lk + (nsec_in_target_time - nsec_in_reftime)

        tsec = real(delta_sec, rk)
    end function seconds_since_datetime

    !-------------------------------------------
    ! Convert seconds since reference date to date/time
    !-------------------------------------------
    pure subroutine seconds_to_datetime(cal, u, tsec, year, month, day, hour, minute, second, doy)
        type(CFCalendar), intent(in) :: cal              ! Calendar
        type(CFUnits),    intent(in) :: u                ! Time units
        real(rk),         intent(in) :: tsec             ! Number of seconds since reference date
        integer,          intent(out):: year, month, day, hour, minute, second, doy

        integer(lk) :: total_s, days_since_ref, remaining_seconds
        integer(lk) :: numdays_in_ref, numdays_in_target_date

        total_s = nint(tsec, lk)
        days_since_ref = total_s / 86400_lk              ! Number of whole days since reference date (Integer)
        remaining_seconds = mod(total_s, 86400_lk)       ! Remaining number of seconds in the last day
        if (remaining_seconds < 0_lk) then
            remaining_seconds = remaining_seconds + 86400_lk
            days_since_ref = days_since_ref - 1_lk
        end if

        numdays_in_ref  = date_to_number_of_days(cal%kind, u%reference_year, u%reference_month, u%reference_day)
        numdays_in_target_date = numdays_in_ref + days_since_ref

        ! Obtains the date for the number of seconds
        call numdays_to_date(cal%kind, numdays_in_target_date, year, month, day)
        hour = int(remaining_seconds / 3600_lk)                  ! Number of hours
        minute = int(mod(remaining_seconds, 3600_lk) / 60_lk)    ! Number of minutes
        second = int(mod(remaining_seconds, 60_lk))              ! Number of seconds
        doy = get_day_of_year(cal%kind, year, month, day)        ! Day of year
    end subroutine seconds_to_datetime

    !----------------------------------------------------------
    ! Number of Days in a given year depending on the Calendar
    !--------------------------------------------------------
    pure integer function get_days_in_year(cal, year) result(nd)
        type(CFCalendar), intent(in) :: cal
        integer, intent(in) :: year
        select case (cal%kind)
        case (cal_gregorian, cal_proleptic)
            ! 366 if it is leap (Gregorian) otherwise 365
            nd = merge(366, 365, is_leap_gregorian(year))
        case (cal_noleap)
            nd = 365
        case (cal_all_leap)
            nd = 366
        case (cal_360_day)
            nd = 360
        case default
            nd = 365
        end select
    end function get_days_in_year



!=======================
! Internals
!======================

    ! Split "time-unit since referemce date"
    subroutine split_units_since(units_attr, time_unit, reference_date_str, ok)
        character(*), intent(in)  :: units_attr
        character(len=:), allocatable, intent(out) :: time_unit, reference_date_str
        logical, intent(out) :: ok
        character(len=:), allocatable :: s
        integer :: p

        s  = trim(units_attr)
        p  = index(to_lower(s), 'since', back=.false.)
        if (p <= 1) then
            ok = .false.; return
        end if
        time_unit  = trim(adjustl(s(1:p-1)))
        reference_date_str = trim(adjustl(s(p+5:)))   ! 5 chars in "since"
        ok = (len_trim(time_unit) > 0 .and. len_trim(reference_date_str) > 0)
    end subroutine split_units_since

    ! Verifies the time-units and assigns number of seconds per time-unit (e.g. 3600 seconds per hour)
    subroutine verify_time_unit(time_unit, timeunit_to_seconds, ok)
        character(*), intent(in)  :: time_unit
        real(rk),    intent(out)  :: timeunit_to_seconds
        logical,     intent(out)  :: ok
        character(len=:), allocatable :: u
        u = to_lower(trim(time_unit))
        select case (u)
        case ('second','seconds','sec','secs','s')
            timeunit_to_seconds = 1.0_rk; ok = .true.
        case ('minute','minutes','min','mins')
            timeunit_to_seconds = 60.0_rk; ok = .true.
        case ('hour','hours','hr','hrs','h')
            timeunit_to_seconds = 3600.0_rk; ok = .true.
        case ('day','days','d')
            timeunit_to_seconds = 86400.0_rk; ok = .true.
        case ('week','weeks')
            timeunit_to_seconds = 7.0_rk*86400.0_rk; ok = .true.
        case ('millisecond','milliseconds','msec','msecs','ms')
            timeunit_to_seconds = 1.0e-3_rk; ok = .true.
        case ('microsecond','microseconds','usec','usecs','us')
            timeunit_to_seconds = 1.0e-6_rk; ok = .true.
        case default
            timeunit_to_seconds = 0.0_rk; ok = .false.
        end select
    end subroutine verify_time_unit

    ! verifies calendar string
    subroutine verify_calendar_str(cal_str, cal_kind, ok)
        character(*), intent(in) :: cal_str
        integer, intent(out) :: cal_kind
        logical, intent(out) :: ok
        select case (to_lower(trim(cal_str)))
        case ('gregorian','standard','')
            cal_kind = cal_gregorian; ok = .true.
        case ('proleptic_gregorian')
            cal_kind = cal_proleptic; ok = .true.
        case ('noleap','365_day')
            cal_kind = cal_noleap; ok = .true.
        case ('all_leap','366_day')
            cal_kind = cal_all_leap; ok = .true.
        case ('360_day')
            cal_kind = cal_360_day; ok = .true.
        case default
            ok = .false.
        end select
    end subroutine verify_calendar_str  
    

    !==============================================
    ! Conversions between dates and number of days
    !==============================================

    ! Calculates the number of days in a given date since the Common Era started
    ! This calculation depends on the kiven calendar. 
    pure integer(lk) function date_to_number_of_days(cal_kind, y, m, d) result(num_days)
        integer, intent(in) :: cal_kind, y, m, d
        integer :: mm
        integer(lk) :: days_before_year

        select case (cal_kind)
        case (cal_360_day)
            num_days = int((y-1), lk)*360_lk + int((m-1)*30 + (d-1), lk)
            return
        case (cal_noleap)
            days_before_year = int(y-1, lk) * 365_lk
            mm = max(1,min(12,m))
            num_days = days_before_year + int(cum_nonleap(mm),lk) + int(d-1, lk)
            return
        case (cal_all_leap)
            days_before_year = int(y-1, lk) * 366_lk
            mm = max(1,min(12,m))
            num_days = days_before_year + int(cum_leap(mm),lk) + int(d-1, lk)
            return
        case (cal_gregorian, cal_proleptic)
            days_before_year = int(365_lk*(y-1) + (y-1)/4 - (y-1)/100 + (y-1)/400, lk)
            mm = max(1,min(12,m))
            if (is_leap_gregorian(y)) then
                num_days = days_before_year + int(cum_leap(mm),lk) + int(d-1, lk)
            else
                num_days = days_before_year + int(cum_nonleap(mm),lk) + int(d-1, lk)
            end if
            return
        case default
            ! Fallback to gregorian
            days_before_year = int(365_lk*(y-1) + (y-1)/4 - (y-1)/100 + (y-1)/400, lk)
            mm = max(1,min(12,m))
            if (is_leap_gregorian(y)) then
                num_days = days_before_year + int(cum_leap(mm),lk) + int(d-1, lk)
            else
                num_days = days_before_year + int(cum_nonleap(mm),lk) + int(d-1, lk)
            end if
        end select
    end function date_to_number_of_days

    ! Convert an integer number of days (num_days) since the Common Era started to a date (Y,M,D)
    ! Assigns year (y), month (m) and day (d)
    pure subroutine numdays_to_date(cal_kind, num_days, y, m, d)        
        integer, intent(in) :: cal_kind
        integer(lk), intent(in) :: num_days
        integer, intent(out) :: y, m, d
        integer(lk) :: days, diy

        select case (cal_kind)
        case (cal_360_day)
            y = int(num_days / 360_lk) + 1
            diy = mod(num_days, 360_lk)
            m = int(diy / 30_lk) + 1
            d = int(mod(diy, 30_lk)) + 1
            return
        case (cal_noleap)
            y = int(num_days / 365_lk) + 1
            diy = mod(num_days, 365_lk)
            call month_day_from_doy_noleap(int(diy), m, d)
            return
        case (cal_all_leap)
            y = int(num_days / 366_lk) + 1
            diy = mod(num_days, 366_lk)
            call month_day_from_doy_leap(int(diy), m, d)
            return
        case (cal_gregorian, cal_proleptic)
            ! Invert days_before_year formula via search (fast enough for realistic ranges)
            y = 1
            ! First, rough guess:
            y = max(1, int( real(num_days, rk) / 365.2425_rk ) + 1)
            do
                days = 365_lk*(y-1) + (y-1)/4 - (y-1)/100 + (y-1)/400
                if (days > num_days) then
                y = y - 1
                exit
                end if
                if (365_lk*y + y/4 - y/100 + y/400 > num_days) exit
                y = y + 1
            end do
            diy = num_days - (365_lk*(y-1) + (y-1)/4 - (y-1)/100 + (y-1)/400)
            if (is_leap_gregorian(y)) then
                call month_day_from_doy_generic(int(diy), cum_leap, m, d)
            else
                call month_day_from_doy_generic(int(diy), cum_nonleap, m, d)
            end if
            return
        case default
            y = int(num_days / 365_lk) + 1
            diy = mod(num_days, 365_lk)
            call month_day_from_doy_noleap(int(diy), m, d)
        end select
    end subroutine numdays_to_date

    ! Returns the day of year for a given date
    pure integer function get_day_of_year(cal_kind, y, m, d) result(doy)
        integer, intent(in) :: cal_kind, y, m, d

        select case (cal_kind)
        case (cal_360_day)
            doy = (m-1)*30 + d
        case (cal_noleap)
            doy = cum_nonleap(m) + d
        case (cal_all_leap)
            doy = cum_leap(m) + d
        case (cal_gregorian, cal_proleptic)
            if (is_leap_gregorian(y)) then
                doy = cum_leap(m) + d
            else
                doy = cum_nonleap(m) + d
            end if
        case default
            doy = cum_nonleap(m) + d
        end select
    end function get_day_of_year

    ! Gets month and day from day of year in a non-leap year
    pure subroutine month_day_from_doy_noleap(diy, m, d)
        integer, intent(in) :: diy ! 0-based days into year
        integer, intent(out):: m, d
        call month_day_from_doy_generic(diy, cum_nonleap, m, d)
    end subroutine month_day_from_doy_noleap

    ! Gets month and day from day of year in a leap year
    pure subroutine month_day_from_doy_leap(diy, m, d)
        integer, intent(in) :: diy
        integer, intent(out):: m, d
        call month_day_from_doy_generic(diy, cum_leap, m, d)
    end subroutine month_day_from_doy_leap

    pure subroutine month_day_from_doy_generic(diy, cum, m, d)
        integer, intent(in) :: diy
        integer, dimension(12), intent(in) :: cum
        integer, intent(out):: m, d
        integer :: mm
        do mm=1,12
            if (diy < cum(mm) + (cum(min(12,mm+1)) - cum(mm))) then
                m = mm
                d = diy - cum(mm) + 1
                return
            end if
        end do
        m = 12
        d = diy - cum(12) + 1
    end subroutine month_day_from_doy_generic

end module cf_time_utils
