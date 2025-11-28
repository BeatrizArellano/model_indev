module time_types
  use precision_types, only: rk, ik     ! Importing real64 and lk
  implicit none
  private

 type, public :: DateTime
    integer(ik) :: year   = 0
    integer(ik) :: month  = 0
    integer(ik) :: day    = 0
    integer(ik) :: hour   = 0
    integer(ik) :: minute = 0
    integer(ik) :: second = 0
    logical     :: has_time = .false.
  end type DateTime

  !--------------------------------------------------------------------
  ! Calendar identifiers, names according to CF Metadata convention
  !--------------------------------------------------------------------
  integer, parameter, public :: cal_gregorian        = 1  ! "gregorian" or "standard"
  integer, parameter, public :: cal_proleptic        = 2  ! "proleptic_gregorian"
  integer, parameter, public :: cal_noleap           = 3  ! "noleap" or "365_day"
  integer, parameter, public :: cal_all_leap         = 4  ! "all_leap" or "366_day"
  integer, parameter, public :: cal_360_day          = 5  ! "360_day"
  integer, parameter, public :: cal_standard         = 6  !  equal to "gregorian"

  type, public :: CFCalendar
    ! kind : one of the calendar IDs
     integer :: kind = cal_gregorian   ! default
  contains
     procedure, pass(self) :: name => cfcalendar_name   ! method
  end type CFCalendar

  type, public :: CFUnits
    ! Holds the CF units since the reference point in time (date)
     real(rk) :: timeunit_to_seconds = 0.0_rk               !conversion factor from the declared CF time unit to seconds (e.g. 86400 for "days since")
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

  public :: calendar_compatible

contains

    pure function cfcalendar_name(self) result(name)
      class(CFCalendar), intent(in) :: self
      character(len=:), allocatable :: name
      select case (self%kind)
      case (cal_gregorian); name = 'gregorian'
      case (cal_proleptic); name = 'proleptic_gregorian'
      case (cal_noleap);    name = 'noleap'            ! aka 365_day
      case (cal_all_leap);  name = 'all_leap'          ! aka 366_day
      case (cal_360_day);   name = '360_day'
      case default;         name = 'unknown'
      end select
    end function cfcalendar_name
   

    logical function calendar_compatible(cala, calb) result(ok)
      implicit none
      integer, intent(in) :: cala, calb

      logical :: a_greg, b_greg
      a_greg = (cala == cal_gregorian) .or. (cala == cal_proleptic) .or. (cala == cal_standard)
      b_greg = (calb == cal_gregorian) .or. (calb == cal_proleptic) .or. (calb == cal_standard)

      if (a_greg .and. b_greg) then
          ok = .true.
      else
          ok = (cala == calb)
      end if
    end function calendar_compatible
  
end module time_types
