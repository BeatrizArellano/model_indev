module calendar_types
  implicit none
  private

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
  
end module calendar_types
