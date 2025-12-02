module test_load_forcing
  use precision_types,  only: rk, lk
  use read_config_yaml, only: ConfigParams
  use geo_utils,        only: LocationInfo
  use time_types,       only: DateTime, CFCalendar
  use load_forcing,     only: ForcingState, ForcingYearData, scan_and_init_forcing, load_year_data, print_forcing_summary
  use tidal_readers, only: read_tidal_parameters
  use tidal,            only: TidalSet, TidalConstituent
  use validation_utils, only: validate_input_dates, validate_location_input
  
  implicit none

  type(TidalSet)     :: Tides
  type(TidalConstituent)     :: m2, s2, k1, o1, n2
  type(ConfigParams)    :: main_cfg, physics_cfg
  type(LocationInfo)    :: location
  type(DateTime)        :: start_datetime, end_datetime
  type(CFCalendar)      :: calendar
  type(ForcingState)    :: FS
  type(ForcingYearData) :: surf

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
                trim(Tides%c(i)%name), Tides%c(i)%smaj, Tides%c(i)%smin, Tides%c(i)%theta, Tides%c(i)%phase
        end do

        do y = start_datetime%year, end_datetime%year
            k = y - FS%sim_y_start + 1

            call load_year_data(FS, k, surf, ok, errmsg)
            if (.not. ok) stop trim(errmsg)

            nt = size(surf%air_temp%data)
            if (size(surf%slp%data) /= nt .or. size(surf%rel_hum%data) /= nt .or. &
                size(surf%short_rad%data) /= nt .or. size(surf%long_rad%data) /= nt .or. &
                size(surf%wind_u10%data) /= nt .or. size(surf%wind_v10%data) /= nt) then
            stop 'Forcing variables have inconsistent lengths'
            end if

            write(fname,'(A,I0,A)') 'forcing_', y, '.dat'
            open(newunit=iu, file=trim(fname), status='replace', action='write', iostat=ios)
            if (ios /= 0) stop 'Failed to open output file'

            write(iu,'(A)') 'surf_air_temp sl_pressure relative_humidity shortwave_radiation longwave_radiation wind_speed wind_v10ection co2_air'
            do i = 1, nt
            write(iu,'(8(1X,ES16.8))') surf%air_temp%data(i), surf%slp%data(i),       &
                                        surf%rel_hum%data(i),  surf%short_rad%data(i), &
                                        surf%long_rad%data(i), surf%wind_u10%data(i),  &
                                        surf%wind_v10%data(i), surf%co2_air%const_value
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
        write(*,'(A,1X,ES12.5)') 'wind_u10:', surf%wind_u10%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'wind_v10:', surf%wind_v10%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'shortwave:', surf%short_rad%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'longwave:', surf%long_rad%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'rel_hum:', surf%rel_hum%value_at_step(time)
        write(*,'(A,1X,ES12.5)') 'co2:', surf%co2_air%value_at_step(time)
    end subroutine test_load_data
end module test_load_forcing


 !integer(lk) :: i, t_step_sec, dt_win

        !!----------- For testing--------------------------------------------
        !integer, parameter :: SECS_PER_DAY = 86400
        !integer, parameter :: n_tests = 6
        !integer, parameter :: test_days(n_tests) = [1, 7, 157, 360, 361, 500]
        !logical            :: printed(n_tests)
        !integer            :: k, day_num
        !printed = .false.
        !!---------------------------------------------------------------------


        !if (.not. is_main_initialized) error stop 'run_shelfseas: shelfseas not initialised.'
        !do i = 1_lk, n_steps
        !    ! Model time at the start of this step (seconds since start_datetime)
        !    t_step_sec = (i - 1_lk) * dt
        !    call ForcMan%tick(t_step_sec)                    ! Time-manager to load yearly forcing data on time
        !    call ForcMan%sample(t_step_sec, ForcSnp)         ! get forcing snapshot for the current model time

        !    !-------------- TEST --------------------------------------------------------
        !    ! Determine (1-based) day number for current step
        !    !day_num = int(t_step_sec / SECS_PER_DAY, kind=lk) + 1
        !    
        !    ! If this is one of the target days and we haven't printed yet, dump the snapshot
        !    do k = 1, n_tests
        !        day_num = ((test_days(k) -1)* SECS_PER_DAY) + (3600 * 11)
        !        !if (.not. printed(k) .and. day_num == test_days(k)) then
        !        if (.not. printed(k) .and. t_step_sec >= day_num) then
        !            write(*,'(A,I0)') '--- Forcing snapshot at day ', test_days(k)
        !            write(*,'(A,F12.5)') '  air_temp       = ', ForcSnp%air_temp
        !            write(*,'(A,F12.5)') '  slp            = ', ForcSnp%slp
        !            write(*,'(A,F12.5)') '  rel_hum        = ', ForcSnp%rel_hum
        !            write(*,'(A,F12.5)') '  short_rad      = ', ForcSnp%short_rad
        !            write(*,'(A,F12.5)') '  long_rad       = ', ForcSnp%long_rad
        !            write(*,'(A,F12.5)') '  wind_u10       = ', ForcSnp%wind_u10
        !            write(*,'(A,F12.5)') '  wind_v10       = ', ForcSnp%wind_v10
        !            write(*,'(A,F12.5)') '  co2_air        = ', ForcSnp%co2_air
        !            printed(k) = .true.
        !        end if
        !    end do
            !--------------------------------------------------------------------------

            ! Check Forcing: ensure active data covers [t_step_sec, t_step_sec + dt_win]            

            ! Physics step (handles its own subcycling/implicit solves as needed)
            ! call physics_step(t_step_sec, dt_win, wgrid, ...)

            ! Transport step (only if enabled; also subcycles internally)
            ! call transport_step(t_step_sec, dt_win, wgrid, ...)

            ! Output/diagnostics if due
            ! call output_maybe_write(t_step_sec, dt_win, ...)