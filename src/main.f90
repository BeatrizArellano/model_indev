program shelfmodel_main
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_utils,       only: DateTime, datetime_to_str
  use calendar_types,   only: CFCalendar
  use validation_utils, only: validate_input_dates, validate_location_input, print_header
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
    ! Print header for simulation
    call print_header(location,start_datetime,end_datetime) 


    call physics_init(location,start_datetime, end_datetime, calendar)

    call physics_end()

    call main_cfg%clear()    

end program shelfmodel_main




