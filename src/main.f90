program shelfmodel_main
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_utils,       only: DateTime
  use calendar_types,   only: CFCalendar
  use validation_utils, only: validate_input_dates, validate_location_input
  use physics_driver

  implicit none

    type(ConfigParams) :: main_cfg
    type(LocationInfo) :: location
    type(DateTime)     :: start_datetime, end_datetime
    type(CFCalendar)   :: calendar

    real(kind=8) :: dt, t
    
    call main_cfg%init()    
    call main_cfg%load_yaml_content('main.yaml')

    ! Read and validate location parameters
    call validate_location_input(main_cfg, location)
    ! Read and validate simulation dates and calendar
    call validate_input_dates(main_cfg, start_datetime, end_datetime, calendar)

    write(*,*) 'Name: ', location%name
    write(*,*) 'latitude=', location%lat, ' longitude=', location%lon, ' depth=', location%depth

    
    write(*,*) 'Calendar: ', trim(adjustl(calendar%name()))
    if (start_datetime%has_time) then
      write(*,'("start: ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
          start_datetime%year, start_datetime%month, start_datetime%day, &
          start_datetime%hour, start_datetime%minute, start_datetime%second
    else
      write(*,'("start: ",I4.4,"-",I2.2,"-",I2.2)') &
          start_datetime%year, start_datetime%month, start_datetime%day
    end if


    call physics_init(location,start_datetime, end_datetime, calendar)

    call physics_end()

    call main_cfg%clear()    

end program shelfmodel_main




