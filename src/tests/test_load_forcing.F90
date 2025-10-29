module test_load_forcing
  use precision_types,  only: rk, lk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_types,       only: DateTime, CFCalendar
  use load_forcing,     only: ForcingState, ForcingYearData, scan_and_init_forcing, load_year_data, print_forcing_summary
  use tidal_parameters_readers, only: TidalParams, Constituent, read_tidal_parameters
  use validation_utils, only: validate_input_dates, validate_location_input
  
  implicit none

  type(TidalParams)  :: Tides
  type(Constituent)  :: m2, s2, k1, o1, n2
  type(ConfigParams) :: main_cfg, physics_cfg
  type(LocationInfo) :: location
  type(DateTime)     :: start_datetime, end_datetime
  type(CFCalendar)   :: calendar
  type(ForcingState) :: FS
  type(ForcingYearData):: surf

  contains
    subroutine test_load_data()

        integer(lk) :: time
        logical :: ok
        character(len=256) :: errmsg
        integer  :: i, y, k, nt, iu, ios
        character(len=128) :: fname

        call main_cfg%init()
        call main_cfg%load_yaml_content('main.yaml')

        call validate_location_input(main_cfg, location)
        call validate_input_dates(main_cfg, start_datetime, end_datetime, calendar)

        call physics_cfg%init()
        call physics_cfg%load_yaml_content('physics.yaml')

        write(*,*) '--- Forcing scan ---'
        call scan_and_init_forcing(physics_cfg, calendar, location, start_datetime, end_datetime, FS, ok, errmsg)
        if (.not. ok) stop trim(errmsg)
        call print_forcing_summary(FS)

        call read_tidal_parameters(physics_cfg, location%lat, location%lon, Tides)
        write(*,*) 'Tidal constituents loaded:', size(Tides%c)
        do i = 1, size(Tides%c)
            write(*,'(A3,2X,"SEMA=",F8.3,2X,"SEMI=",F8.3,2X,"INC=",F7.2," deg",2X,"PHA=",F7.2," deg")') &
                trim(Tides%c(i)%name), Tides%c(i)%sema, Tides%c(i)%semi, Tides%c(i)%inc_deg, Tides%c(i)%pha_deg
        end do

        do y = start_datetime%year, end_datetime%year
            k = y - FS%sim_y_start + 1

            call load_year_data(FS, k, surf, ok, errmsg)
            if (.not. ok) stop trim(errmsg)

            nt = size(surf%air_temp%data)
            if (size(surf%slp%data) /= nt .or. size(surf%rel_hum%data) /= nt .or. &
                size(surf%short_rad%data) /= nt .or. size(surf%long_rad%data) /= nt .or. &
                size(surf%wind_spd%data) /= nt .or. size(surf%wind_dir%data) /= nt) then
            stop 'Forcing variables have inconsistent lengths'
            end if

            write(fname,'(A,I0,A)') 'forcing_', y, '.dat'
            open(newunit=iu, file=trim(fname), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'Failed to open output file'

            write(iu,'(A)') 'surf_air_temp sl_pressure relative_humidity shortwave_radiation longwave_radiation wind_speed wind_direction co2_air'
            do i = 1, nt
            write(iu,'(8(1X,ES16.8))') surf%air_temp%data(i), surf%slp%data(i),       &
                                        surf%rel_hum%data(i),  surf%short_rad%data(i), &
                                        surf%long_rad%data(i), surf%wind_spd%data(i),  &
                                        surf%wind_dir%data(i), surf%co2_air%const_value
            end do
            close(iu)
        end do

        call load_year_data(FS, 1, surf, ok, errmsg)
        if (.not. ok) stop trim(errmsg)

        time = 30_lk*86400_lk + 10_lk
        write(*,'(A,1X,ES12.5)') 'Tair_:',   surf%air_temp%value_at_step(time)
        time = 350_lk*86400_lk + 86300_lk
        write(*,'(A,1X,ES12.5)') 'Tair_180:', surf%air_temp%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'psl:', surf%slp%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'wind_speed:', surf%wind_spd%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'wind_dir:', surf%wind_dir%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'shortwave:', surf%short_rad%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'longwave:', surf%long_rad%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'rel_hum:', surf%rel_hum%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'co2:', surf%co2_air%value_at_step(time)
    end subroutine test_load_data
end module test_load_forcing