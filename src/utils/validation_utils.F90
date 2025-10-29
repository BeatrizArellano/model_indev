module validation_utils
  use precision_types,  only: rk, ik    ! Importing real64 and int32
  use precision_utils,  only: round_to
  use read_config_yaml, only: ConfigParams
  use time_types,       only: DateTime, CFCalendar, cal_gregorian
  use time_utils,       only: datetime_to_str, &                        
                              validate_calendar, datetime_from_string, is_datetime_before, is_datetime_equal
  use geo_utils,        only: LocationInfo, is_lat_valid, is_lon_valid, is_depth_valid, convert_to_lon180
  use str_utils,        only: realtostr

  implicit none
  private

  public :: validate_input_dates, validate_location_input, print_header

contains

    ! Read and validate dates + calendar from YAML.
    ! ERROR STOPs if:
    !  - calendar code is invalid,
    !  - dates can’t be parsed/validated for that calendar,
    !  - or start >= end.
    subroutine validate_input_dates(user_cfg, start_datetime, end_datetime, calendar)
        type(ConfigParams), intent(in)  :: user_cfg
        type(DateTime),     intent(out) :: start_datetime, end_datetime
        type(CFCalendar),   intent(out) :: calendar

        character(len=:), allocatable :: start_date_str, end_date_str
        integer :: calcode, calcode_val
        logical :: ok
        character(len=256) :: msg
        character(len=:), allocatable :: stopmsg

        ! Read dates and calendar from yaml config files
        start_date_str = user_cfg%get_param_str('time.start', required=.true., trim_value=.true.)
        end_date_str   = user_cfg%get_param_str('time.end', required=.true., trim_value=.true.)
        calcode        = user_cfg%get_param_int('time.calendar', default=0, min=0, max=5)


        ! Validate calendar
        call validate_calendar(calcode, ok)
        if (calcode == 0) then
            calendar%kind = 0
            calcode_val = cal_gregorian              ! Assume Gregorian for validations
        else if (ok) then
            calendar%kind = calcode
            calcode_val   = calcode
        else
            write(msg,'(A,I0,A)') 'Invalid calendar code ', calcode, ' (valid: 0 to 5).'
            error stop trim(msg)
        end if
        
        ! Set CFCalendar
        calendar%kind = calcode

        !Parse & validate dates under that calendar 
        start_datetime = datetime_from_string(trim(start_date_str), ok, msg, calcode_val)
        end_datetime   = datetime_from_string(trim(end_date_str),   ok, msg, calcode_val)

        ! Normalize missing times: start -> 00:00:00, end -> 23:59:59 (inclusive)
        if (.not. start_datetime%has_time) then
            start_datetime%hour   = 0
            start_datetime%minute = 0
            start_datetime%second = 0
            start_datetime%has_time = .false.
        end if

        if (.not. end_datetime%has_time) then
            end_datetime%hour   = 23
            end_datetime%minute = 59
            end_datetime%second = 59
            end_datetime%has_time = .false.
        end if


        ! Validate that start_datetime is before end_datetime
        if (.not. is_datetime_before(start_datetime, end_datetime)) then
            if (is_datetime_equal(start_datetime, end_datetime)) then
                stopmsg = 'Invalid simulation time period: start datetime equals end ('// &
                        trim(start_date_str)//' == '//trim(end_date_str)//').'
            else
                stopmsg = 'Invalid simulation time period: start datetime ('//trim(start_date_str)// &
                        ') is not before end ('//trim(end_date_str)//').'
            end if
            error stop stopmsg
        end if
    end subroutine validate_input_dates


    !-----------------------------------------------------------------------
    ! Read and validate location block:
    !   main.location.name      (optional string)
    !   main.location.latitude  (required, degrees North: -90..90)
    !   main.location.longitude (required, degrees East: -180..180 or 0..360)
    !   main.location.depth     (required, meters > 0)
    !
    ! Returns:
    !   loc : LocationInfo with lon normalized to [-180,180)
    !
    ! ERROR STOPs with a clear message if any field is invalid.
    !-----------------------------------------------------------------------
    subroutine validate_location_input(user_cfg, loc)
        type(ConfigParams), intent(in)  :: user_cfg
        type(LocationInfo), intent(out) :: loc

        character(len=:), allocatable :: name_str
        real(rk) :: lat, lon, depth
        character(len=160) :: msg

        ! Read fields (name can be absent/empty depending on your API)
        name_str = user_cfg%get_param_str('location.name', default='', empty_ok=.true., trim_value=.true.)
        lat      = user_cfg%get_param_num('location.latitude',  required=.true.)
        lon      = user_cfg%get_param_num('location.longitude', required=.true.)
        depth    = user_cfg%get_param_num('location.depth', required=.true.)

        ! Validate latitude
        if (.not. is_lat_valid(lat)) then
        write(msg,'(A,F0.6,A)') 'Invalid latitude: ', lat, ' (expected -90 to 90 decimal degrees).'
        error stop trim(msg)
        end if

        ! Validate longitude
        if (.not. is_lon_valid(lon)) then
        write(msg,'(A,F0.6,A)') 'Invalid longitude: ', lon, ' (expected -180..180 or 0..360 decimal degrees).'
        error stop trim(msg)
        end if

        ! Validate depth
        if (.not. is_depth_valid(depth)) then
        write(msg,'(A,F0.3)') 'Invalid depth (must be > 0 m): ', depth
        error stop trim(msg)
        end if

        ! Assign data to loc and convert longitude to [-180,180)
        loc%lat   = lat
        loc%lon   = convert_to_lon180(lon)
        loc%depth = depth
        loc%name = trim(name_str)
    end subroutine validate_location_input

    subroutine print_header(location, start_datetime, end_datetime)
        use, intrinsic :: iso_fortran_env, only: output_unit
        implicit none
        type(LocationInfo), intent(in) :: location
        type(DateTime),     intent(in) :: start_datetime, end_datetime

        character(len=32) :: slat, slon, sdep

        ! Optional location name (only if non-empty)
        if (len_trim(location%name) > 0) then
            write(output_unit,'(7X,"Location: ",A)') trim(location%name)
        end if
        write(output_unit,'("Latitude: ",F0.3,4X,"Longitude:",F0.1,4X,"Depth: ",F0.1," m")') &
                            location%lat, location%lon, location%depth


        ! Simulation start / end
        write(output_unit,'("start: ",A,"   end: ",A)') trim(datetime_to_str(start_datetime)), &
                                                            trim(datetime_to_str(end_datetime))
    end subroutine print_header


end module validation_utils