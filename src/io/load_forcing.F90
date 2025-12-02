module load_forcing
   use precision_types,     only: rk, lk
   use str_utils,           only: to_lower, inttostr
   use read_config_yaml,    only: ConfigParams
   use geo_utils,           only: LocationInfo
   use time_types,          only: DateTime, CFCalendar, calendar_compatible
   use read_forcing_ncdata, only: ForcingScan, scan_forcing, get_ts_slice_from_point
   use cf_time_utils,       only: seconds_since_datetime_file
   use time_utils,          only: index_at_or_before, index_at_or_after
   use netcdf_io,           only: NcFile, nc_open, nc_close

   implicit none
   private
   public :: ForcingCfg, VarInfo, ForcingState, ForcingYearData, ForcingVarData, scan_and_init_forcing, load_year_data, print_forcing_summary

   integer,     parameter :: in_file=1, in_constant=2, in_compute=3, in_off=4
   integer(lk), parameter :: INF_EDGE = huge(1_lk)

   type :: ForcingCfg
      logical :: repeat_enabled = .false.
      integer :: repeat_year    = -huge(1)
      integer :: cfg_calendar = 0  
   end type ForcingCfg

   type :: VarInfo
      character(:), allocatable :: id                    ! Variable internal name
      integer                   :: input_type            ! file, constant, compute, off
      character(:), allocatable :: name_in_file
      character(:), allocatable :: file_path             ! per-var or global; empty if constant/compute
      integer                   :: calendar    
      real(rk)                  :: first_time, last_time
      integer, allocatable      :: idx_year(:,:)         ! Time indices for when a year starts and ends [start_yi, end_yi]

      ! link to file scan
      type(ForcingScan), pointer :: scan => null()       ! time axis, units, etc.

      ! If the value is constant. 
      real(rk)                  :: const_value = 0.0_rk
   end type VarInfo

   type :: ForcingState
      type(ForcingCfg)               :: cfg
      type(VarInfo),     allocatable :: vars(:)            ! per-variable plan (idx_year etc.)
      character(:),      allocatable :: files(:)           ! unique file paths
      type(ForcingScan), allocatable :: scans(:)           ! one per file
      type(CFCalendar)               :: sim_cal
      integer,           allocatable :: scan_idx_of_var(:) ! map var -> scans(:) index (0 if const/off)
      integer                        :: sim_y_start = 0, sim_y_end = 0
   end type

   type :: ForcingVarData
      character(:), allocatable :: name
      logical                   :: is_const = .false.
      real(rk)                  :: const_value = 0._rk
      integer(lk), allocatable  :: t_axis(:)             ! Time-stamps in model time (seconds since start datetime of simulation)
      integer(lk), allocatable  :: t_edge(:)             ! Edges for time steps: midpoints between consecutive stamps. Each edge pair defines the interval where a value is valid.
      real(rk), allocatable     :: data(:)               ! year slice data (omit if const)
      integer                   :: idx   = 1             ! current index in t_axis/data - Initialises the cursor (we start in time interval 1)
      integer                   :: n     = 0
      integer(lk)               :: t_next = huge(1_lk)   ! The boundary when we must switch to the next time-interval.
      contains
         procedure :: init_cursor
         procedure :: value_at_step
   end type

   type :: ForcingYearData
      type(ForcingVarData) :: air_temp
      type(ForcingVarData) :: slp
      type(ForcingVarData) :: rel_hum
      type(ForcingVarData) :: short_rad
      type(ForcingVarData) :: long_rad

      ! Wind components
      type(ForcingVarData) :: wind_u10
      type(ForcingVarData) :: wind_v10

      ! Freshwater fluxes
      type(ForcingVarData) :: precip     ! precipitation
      type(ForcingVarData) :: evap       ! evaporation
      type(ForcingVarData) :: runoff     ! runoff / river input

      type(ForcingVarData) :: co2_air
   end type


   character(len=24), parameter :: var_names(11) = [ character(len=24) :: &
                                                   'surf_air_temp',       'sl_pressure',        'relative_humidity',  &
                                                   'shortwave_radiation', 'longwave_radiation',                        &
                                                   'wind_u10',            'wind_v10',                                  &
                                                   'precipitation',       'evaporation',          'runoff',           &
                                                   'co2_air' ]



contains


   subroutine scan_and_init_forcing(params, calendar_cfg, location, start_datetime, end_datetime, state, ok, errmsg)
      implicit none
      ! inputs
      type(ConfigParams), intent(in)    :: params
      type(CFCalendar),   intent(in)    :: calendar_cfg
      type(LocationInfo), intent(in)    :: location
      type(DateTime),     intent(in)    :: start_datetime, end_datetime
      ! outputs
      type(ForcingState), intent(out)   :: state
      logical,            intent(out)   :: ok
      character(*),       intent(out)   :: errmsg

      ! locals
      type(VarInfo),       allocatable :: vars(:)
      character(:),        allocatable :: files(:)
      type(ForcingScan),   allocatable :: scans(:)
      type(ForcingCfg)                 :: cfg
      character(:),        allocatable :: global_file
      character(:),        allocatable :: req_names(:)
      character(len=8), dimension(2) :: sal_choices
      character(:), allocatable :: sal_mode
      integer :: idx_p, idx_e
      logical :: ok_inside
      integer :: i, j, nfiles, nreq, maxlen, refk

      ok = .false.; errmsg = ''
      ok_inside = .false.

      ! Read forcing parameters
      call read_forcing_config(params, calendar_cfg, start_datetime, end_datetime, vars, cfg, global_file)

      sal_choices = ['constant','compute ']

      sal_mode = to_lower(trim(params%get_param_str('physics.variables.salinity.mode', &
                                                 choices=sal_choices, trim_value=.true., &
                                                 match_case=.false., default='constant')))

      if (sal_mode == 'constant') then
         ! If salinity is constant, completely ignore all freshwater fluxes.
         do i = 1, size(vars)
            select case (trim(vars(i)%id))
            case ('precipitation','evaporation','runoff')
               if (vars(i)%input_type /= in_off) then
                  vars(i)%input_type = in_off
                  vars(i)%name_in_file = ''
                  vars(i)%file_path    = ''
               end if
            end select
         end do

      else if (sal_mode == 'compute') then
         ! Salinity is prognostic:
         ! - precipitation and evaporation are REQUIRED
         ! - runoff remains optional
         idx_p = find_var_index(vars, 'precipitation')
         idx_e = find_var_index(vars, 'evaporation')

         if (idx_p <= 0 .or. vars(idx_p)%input_type == in_off) then
            errmsg = 'Salinity mode=compute but forcing.precipitation is off or missing. ' // &
                     'Set forcing.precipitation.mode to file or constant, or use salinity.mode=constant.'
            call clear_state(state)
            return
         end if

         if (idx_e <= 0 .or. vars(idx_e)%input_type == in_off) then
            errmsg = 'Salinity mode=compute but forcing.evaporation is off or missing. ' // &
                     'Set forcing.evaporation.mode to file or constant, or use salinity.mode=constant.'
            call clear_state(state)
            return
         end if
      end if

      ! Unique file list
      call build_file_list(vars, files, nfiles)

      ! Build a union list of variable names to validate (simple path)
      if (nfiles > 0) then
         nreq = 0; maxlen = 1
         do i=1, size(vars)
            if (vars(i)%input_type == in_file) then
            nreq   = nreq + 1
            maxlen = max(maxlen, len_trim(vars(i)%name_in_file))
            end if
         end do
         if (nreq == 0) then
            allocate(character(len=1) :: req_names(0))
         else
            allocate(character(len=maxlen) :: req_names(nreq))
            nreq = 0
            do i=1, size(vars)
            if (vars(i)%input_type == in_file) then
               nreq = nreq + 1
               req_names(nreq) = trim(vars(i)%name_in_file)
            end if
            end do
         end if

         
         ! Scan each file once
         call scan_and_attach(files, location%lat, location%lon,                                 &
                              start_datetime%year, start_datetime%month, start_datetime%day,     &
                              start_datetime%hour, start_datetime%minute, start_datetime%second, &
                              end_datetime%year,   end_datetime%month,   end_datetime%day,       &
                              end_datetime%hour,   end_datetime%minute,   end_datetime%second,   &
                              calendar_cfg%name(), req_names, scans, ok_inside, errmsg)
         if (.not. ok_inside) then
            call clear_state(state)
            return
         end if

      else
         allocate(character(len=1) :: files(0))
         allocate(scans(0))
      end if

      ! Attach scans to variables, set calendar, first/last
      call attach_vars_to_scans(vars, files, scans)

      ! Sim year range
      state%sim_y_start = start_datetime%year
      state%sim_y_end   = end_datetime%year

      ! Check calendar compatibility across files
      if (size(scans) > 1) then
         refk = scans(1)%cal%kind
         do i = 2, size(scans)
            if (.not. calendar_compatible(scans(i)%cal%kind, refk)) then
               errmsg = 'Calendar mismatch across forcing files. Preprocess to a common calendar.'
               call clear_state(state)
               return
            end if
         end do
      end if

      ! Decide simulation calendar (derive if cfg=0, else use config)
      if (size(scans) > 0) then
         if (cfg%cfg_calendar == 0) then
            state%sim_cal = scans(1)%cal                       ! derive from forcing
         else
            if (.not. calendar_compatible(cfg%cfg_calendar, scans(1)%cal%kind)) then   
               errmsg = 'Calendar mismatch: forcing uses '//trim(scans(1)%cal%name())// &
                        ', but config is set to '//trim(calendar_cfg%name())//'. ' // &
                        'Set time.calendar: 0 to follow the calendar in the forcing data, or convert your forcing files to match.'
               call clear_state(state)
               return
            end if
            state%sim_cal%kind = cfg%cfg_calendar
         end if
      end if 


      ! Build per-variable year index windows
      do i=1, size(vars)
         if (vars(i)%input_type == in_file) then
            if (.not. associated(vars(i)%scan)) then
               errmsg = 'scan_and_init_forcing: missing scan for '//trim(vars(i)%id)
               call clear_state(state)
               return
            end if
            call build_year_indices(vars(i), state%sim_y_start, state%sim_y_end, cfg%repeat_enabled, cfg%repeat_year)
         end if
      end do

      ! Build var -> scan index map (0 if const/off)
      allocate(state%scan_idx_of_var(size(vars)))
      do i=1, size(vars)
         state%scan_idx_of_var(i) = 0
         if (vars(i)%input_type == in_file) then
            do j=1, size(files)
            if (trim(vars(i)%file_path) == trim(files(j))) then
               state%scan_idx_of_var(i) = j; exit
            end if
            end do
            if (state%scan_idx_of_var(i) == 0) then
               errmsg = 'scan_and_init_forcing: file not found in scans for var '//trim(vars(i)%id)
               call clear_state(state)
               return
            end if
         end if
      end do

      ! Assign to ForcingState
      state%cfg = cfg

      if (allocated(state%vars))  deallocate(state%vars)
      if (allocated(state%files)) deallocate(state%files)
      if (allocated(state%scans)) deallocate(state%scans)

      allocate(state%vars(size(vars)));   state%vars  = vars
      allocate(state%scans(size(scans))); state%scans = scans
      if (size(files) > 0) then
         allocate(character(len=len(files(1))) :: state%files(size(files)))
         state%files = files
      else
         allocate(character(len=1) :: state%files(0))
      end if

      ok = .true.
   end subroutine scan_and_init_forcing

   subroutine read_forcing_config(params, calendar, start_datetime, end_datetime, vars, ctl, global_file)
      type(ConfigParams),  intent(in)      :: params
      type(CFCalendar),    intent(in)         :: calendar
      type(DateTime),      intent(in)         :: start_datetime, end_datetime
      type(VarInfo),       allocatable, intent(out) :: vars(:)
      type(ForcingCfg),        intent(out) :: ctl
      character(:), allocatable,   intent(out) :: global_file

      integer :: i
      character(:), allocatable :: s_mode, s_name, s_fname
      character(len=8), dimension(5) :: mode_choices
      integer :: mode_enum

      mode_choices = ['file    ','constant','compute ','off     ','false   ']
      allocate(vars(size(var_names)))

      ! Global file ('' if missing/null)
      global_file = trim(params%get_param_str('forcing.filename', default=''))

      ! Calendar, we just store kind for now
      ctl%cfg_calendar = calendar%kind

      ! repeat_forcing_year
      if (.not. params%is_disabled('forcing.repeat_forcing_year')) then
         ctl%repeat_year    = params%get_param_int('forcing.repeat_forcing_year')
         ctl%repeat_enabled = .true.
         if (ctl%repeat_year < start_datetime%year .or. ctl%repeat_year > end_datetime%year) then
            stop 'Physics forcing data: repeat_forcing_year: '//inttostr(ctl%repeat_year)// ' is outside simulation time interval.'
         end if
      else
         ctl%repeat_enabled = .false.
         ctl%repeat_year    = -huge(1)
      end if

      do i=1, size(var_names)
         vars(i)%id = var_names(i)

         s_mode = params%get_param_str('forcing.'//trim(adjustl(var_names(i)))//'.mode', &
                                       required=.true., choices=mode_choices, &
                                       trim_value=.true., match_case=.false.)

         select case (trim(s_mode))
         case ('file');     mode_enum = in_file
         case ('constant'); mode_enum = in_constant
         case ('compute');  mode_enum = in_compute
         case ('off','false'); mode_enum = in_off
         end select
         vars(i)%input_type = mode_enum

         select case (mode_enum)
         case (in_file)
            s_name  = params%get_param_str('forcing.'//trim(adjustl(var_names(i)))//'.name',     required=.true.)
            s_fname = params%get_param_str('forcing.'//var_names(i)//'.filename', default='')
            vars(i)%name_in_file = trim(s_name)
            if (len_trim(s_fname) > 0) then
            vars(i)%file_path = trim(s_fname)
            else
            if (len_trim(global_file) == 0) stop 'Var '//trim(var_names(i))// ': mode=file but no filename (per-var or global).'
            vars(i)%file_path = global_file
            end if

         case (in_constant)
            vars(i)%const_value = params%get_param_num('forcing.'//trim(adjustl(var_names(i)))//'.constant', finite=.true.)
            vars(i)%name_in_file = ''     ! not used
            vars(i)%file_path    = ''     ! not used

         case (in_compute)
            stop 'Var '//trim(var_names(i))//': mode=compute not implemented yet.'  ! placeholder
            ! (when enabled later, you may require a name if computed var maps to an output id)

         case (in_off)
            vars(i)%name_in_file = ''
            vars(i)%file_path    = ''
         end select
      end do
   end subroutine read_forcing_config

   subroutine load_year_data(forcmd, year_k, Y, ok, errmsg)
      type(ForcingState),    intent(in)    :: forcmd
      integer,               intent(in)    :: year_k
      type(ForcingYearData), intent(inout) :: Y
      logical,               intent(out)   :: ok
      character(*),          intent(out)   :: errmsg

      type(NcFile) :: db
      integer :: i, j, i0, i1, nt

      ok = .false.; errmsg = ''

      ! -------- Variables from files (open each file once) --------
      do j = 1, size(forcmd%files)
         call nc_open(db, forcmd%files(j))

         do i = 1, size(forcmd%vars)
            if (forcmd%vars(i)%input_type /= in_file)                cycle
            if (forcmd%scan_idx_of_var(i) /= j)                      cycle
            if (.not. allocated(forcmd%vars(i)%idx_year)) then
               errmsg = 'load_year_data: idx_year not built for '//trim(forcmd%vars(i)%id); call nc_close(db); return
            end if
            if (year_k < 1 .or. year_k > size(forcmd%vars(i)%idx_year,2)) then
               errmsg = 'load_year_data: year index out of range for '//trim(forcmd%vars(i)%id); call nc_close(db); return
            end if

            i0 = forcmd%vars(i)%idx_year(1, year_k)
            i1 = forcmd%vars(i)%idx_year(2, year_k)
            nt = max(0, i1 - i0 + 1)
            if (nt <= 0) then
               errmsg = 'load_year_data: empty slice for '//trim(forcmd%vars(i)%id)//' in file '//trim(forcmd%files(j))
               call nc_close(db); return
            end if

            select case (trim(forcmd%vars(i)%id))
            case ('surf_air_temp')
               call load_var_series( Y%air_temp, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%air_temp%name = forcmd%vars(i)%id

            case ('sl_pressure')
               call load_var_series( Y%slp, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%slp%name = forcmd%vars(i)%id

            case ('relative_humidity')
               call load_var_series( Y%rel_hum, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%rel_hum%name = forcmd%vars(i)%id

            case ('shortwave_radiation')
               call load_var_series( Y%short_rad, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%short_rad%name = forcmd%vars(i)%id

            case ('longwave_radiation')
               call load_var_series( Y%long_rad, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%long_rad%name = forcmd%vars(i)%id

            case ('wind_u10')
               call load_var_series( Y%wind_u10, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%wind_u10%name = forcmd%vars(i)%id

            case ('wind_v10')
               call load_var_series( Y%wind_v10, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%wind_v10%name = forcmd%vars(i)%id

            case ('precipitation')
               call load_var_series( Y%precip, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%precip%name = forcmd%vars(i)%id

            case ('evaporation')
               call load_var_series( Y%evap, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%evap%name = forcmd%vars(i)%id

            case ('runoff')
               call load_var_series( Y%runoff, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%runoff%name = forcmd%vars(i)%id


            case ('co2_air')
               call load_var_series( Y%co2_air, forcmd%scans(j), db, forcmd%vars(i)%name_in_file, &
                                    forcmd%scans(j)%time_name, i0, i1, forcmd%scans(j)%has_latlon, &
                                    forcmd%scans(j)%yi, forcmd%scans(j)%xi )
               Y%co2_air%name = forcmd%vars(i)%id

            case default
               errmsg = 'load_year_data: unknown variable id='//trim(forcmd%vars(i)%id)
               call nc_close(db); return
            end select
         end do

         call nc_close(db)
      end do

      ! -------- CONSTANTS --------
      do i = 1, size(forcmd%vars)
         if (forcmd%vars(i)%input_type == in_constant) then
            select case (trim(forcmd%vars(i)%id))
            case ('surf_air_temp')
               Y%air_temp%is_const    = .true.
               Y%air_temp%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%air_temp%t_axis)) deallocate(Y%air_temp%t_axis)
               if (allocated(Y%air_temp%data)) deallocate(Y%air_temp%data)
               Y%air_temp%idx = 1; Y%air_temp%n = 0; Y%air_temp%t_next = INF_EDGE

            case ('sl_pressure')
               Y%slp%is_const    = .true.;  Y%slp%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%slp%t_axis))      deallocate(Y%slp%t_axis)
               if (allocated(Y%slp%data)) deallocate(Y%slp%data)
               Y%slp%idx = 1; Y%slp%n = 0; Y%slp%t_next = INF_EDGE

            case ('relative_humidity')
               Y%rel_hum%is_const = .true.; Y%rel_hum%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%rel_hum%t_axis))  deallocate(Y%rel_hum%t_axis)
               if (allocated(Y%rel_hum%data)) deallocate(Y%rel_hum%data)
               Y%rel_hum%idx = 1; Y%rel_hum%n = 0; Y%rel_hum%t_next = INF_EDGE

            case ('shortwave_radiation')
               Y%short_rad%is_const = .true.; Y%short_rad%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%short_rad%t_axis))deallocate(Y%short_rad%t_axis)
               if (allocated(Y%short_rad%data)) deallocate(Y%short_rad%data)
               Y%short_rad%idx = 1; Y%short_rad%n = 0; Y%short_rad%t_next = INF_EDGE

            case ('longwave_radiation')
               Y%long_rad%is_const = .true.; Y%long_rad%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%long_rad%t_axis)) deallocate(Y%long_rad%t_axis)
               if (allocated(Y%long_rad%data)) deallocate(Y%long_rad%data)
               Y%long_rad%idx = 1; Y%long_rad%n = 0; Y%long_rad%t_next = INF_EDGE

            case ('wind_u10')
               Y%wind_u10%is_const    = .true.
               Y%wind_u10%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%wind_u10%t_axis)) deallocate(Y%wind_u10%t_axis)
               if (allocated(Y%wind_u10%data))   deallocate(Y%wind_u10%data)
               Y%wind_u10%idx = 1; Y%wind_u10%n = 0; Y%wind_u10%t_next = INF_EDGE

            case ('wind_v10')
               Y%wind_v10%is_const    = .true.
               Y%wind_v10%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%wind_v10%t_axis)) deallocate(Y%wind_v10%t_axis)
               if (allocated(Y%wind_v10%data))   deallocate(Y%wind_v10%data)
               Y%wind_v10%idx = 1; Y%wind_v10%n = 0; Y%wind_v10%t_next = INF_EDGE

            case ('precipitation')
               Y%precip%is_const    = .true.
               Y%precip%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%precip%t_axis)) deallocate(Y%precip%t_axis)
               if (allocated(Y%precip%data))   deallocate(Y%precip%data)
               Y%precip%idx = 1; Y%precip%n = 0; Y%precip%t_next = INF_EDGE

            case ('evaporation')
               Y%evap%is_const    = .true.
               Y%evap%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%evap%t_axis)) deallocate(Y%evap%t_axis)
               if (allocated(Y%evap%data))   deallocate(Y%evap%data)
               Y%evap%idx = 1; Y%evap%n = 0; Y%evap%t_next = INF_EDGE

            case ('runoff')
               Y%runoff%is_const    = .true.
               Y%runoff%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%runoff%t_axis)) deallocate(Y%runoff%t_axis)
               if (allocated(Y%runoff%data))   deallocate(Y%runoff%data)
               Y%runoff%idx = 1; Y%runoff%n = 0; Y%runoff%t_next = INF_EDGE

            case ('co2_air')
               Y%co2_air%is_const = .true.; Y%co2_air%const_value = forcmd%vars(i)%const_value
               if (allocated(Y%co2_air%t_axis))  deallocate(Y%co2_air%t_axis)
               if (allocated(Y%co2_air%data)) deallocate(Y%co2_air%data)
               Y%co2_air%idx = 1; Y%co2_air%n = 0; Y%co2_air%t_next = INF_EDGE
            case default
            errmsg = 'load_year_data: unknown constant id='//trim(forcmd%vars(i)%id); return
            end select
         end if
      end do

      ! 
      do i = 1, size(forcmd%vars)
         if (forcmd%vars(i)%input_type == in_off) then
            select case (trim(forcmd%vars(i)%id))
            case ('precipitation')
               Y%precip%is_const    = .true.
               Y%precip%const_value = 0._rk
            case ('evaporation')
               Y%evap%is_const      = .true.
               Y%evap%const_value   = 0._rk
            case ('runoff')
               Y%runoff%is_const    = .true.
               Y%runoff%const_value = 0._rk
            end select
         end if
      end do


      ! Initialize cursors
      call Y%air_temp%init_cursor();  call Y%slp%init_cursor();      call Y%rel_hum%init_cursor()
      call Y%short_rad%init_cursor(); call Y%long_rad%init_cursor()
      call Y%wind_u10%init_cursor();  call Y%wind_v10%init_cursor()
      call Y%precip%init_cursor(); call Y%evap%init_cursor(); call Y%runoff%init_cursor()
      call Y%co2_air%init_cursor()

      ok = .true.
   end subroutine load_year_data


   subroutine init_cursor(self)
      class(ForcingVarData), intent(inout) :: self

      if (self%is_const) then
         self%n = 0; self%idx = 1; self%t_next = INF_EDGE
         return
      end if

      if (allocated(self%t_axis)) then
         self%n = size(self%t_axis)
      else
         self%n = 0
      end if

      self%idx = 1
      if (self%n <= 1) then
         self%t_next = INF_EDGE
      else
         self%t_next = self%t_edge(2)   ! first boundary
      end if
   end subroutine


   
   ! Gets the value at a given t (model time in seconds) for a forcing variable
   real(rk) function value_at_step(self, model_time) result(v)
      class(ForcingVarData), intent(inout) :: self
      integer(lk),           intent(in)    :: model_time

      integer :: nloc

      ! If it has a constant value, return that
      if (self%is_const) then
         v = self%const_value
         return
      end if

      nloc = self%n
      if (.not. self%is_const) then
         if (self%n <= 0) error stop 'value_at_step: empty series for '//trim(self%name)
      end if

      ! Advance when crossing the edge of the current time-step (interval boundary)
      ! If model_time is still before t_next, do nothing
      ! Advancing in a while loop "just in case" the model time-step jumped across more than 1 forcing interval
      do while (self%idx < nloc .and. model_time >= self%t_edge(self%idx+1))
         ! If t has reached or passed t_next, advance. 
         self%idx = self%idx + 1
         if (self%idx < nloc) then
            self%t_next = self%t_edge(self%idx+1)   ! Update t_next to the next edge
         else
            self%t_next = INF_EDGE  ! last interval extends “to infinity” on the right
         end if
      end do

      v = self%data(self%idx)
   end function value_at_step


   subroutine print_forcing_summary(state)
      use, intrinsic :: iso_fortran_env, only: output_unit
      type(ForcingState), intent(in) :: state

      integer :: i
      character(:), allocatable :: list_file, list_const
      logical :: any_file, any_const

      list_file  = ''   ! auto-alloc; we'll concatenate names
      list_const = ''
      any_file   = .false.
      any_const  = .false.

      do i = 1, size(state%vars)
         select case (state%vars(i)%input_type)
         case (in_file)
            if (any_file) then
               list_file = trim(list_file)//', '//trim(state%vars(i)%id)
            else
               list_file = trim(state%vars(i)%id)
               any_file = .true.
            end if
         case (in_constant)
            if (any_const) then
               list_const = trim(list_const)//', '//trim(state%vars(i)%id)
            else
               list_const = trim(state%vars(i)%id)
               any_const = .true.
            end if
         end select
      end do

      if (any_file)   write(*,'(A)') '  ✓ Forcing files scanned - all checks passed'
      if (any_file)   write(output_unit,'(2X,"Variables from file: ",A)')     trim(list_file)
      if (any_const)  write(output_unit,'(2X,"Constant variables: ",A)') trim(list_const) 
      write(*,'(2X,A,A)') 'Calendar: ', trim(state%sim_cal%name())
   end subroutine

   !---------------- HELPERS --------------------------------------------


   ! Builds a unique list of filepaths
   subroutine build_file_list(vars, file_paths, nfiles)
      type(VarInfo),             intent(in)  :: vars(:)
      character(:), allocatable, intent(out) :: file_paths(:)
      integer,                   intent(out) :: nfiles

      integer :: i, j, maxlen
      logical :: found
      character(:), allocatable :: tmp(:)
      character(:), allocatable :: fname

      ! Max length across paths to allocate character(len=...) arrays
      maxlen = 1
      do i=1, size(vars)
         if (vars(i)%input_type == in_file) then
            maxlen = max(maxlen, len_trim(vars(i)%file_path))
         end if
      end do

      allocate(character(len=maxlen) :: tmp(size(vars)))
      tmp = '' ; nfiles = 0

      do i=1, size(vars)
         if (vars(i)%input_type /= in_file) cycle
         fname = trim(vars(i)%file_path)
         if (len_trim(fname) == 0) cycle

         found = .false.
         do j=1, nfiles
            if (trim(tmp(j)) == fname) then
            found = .true.; exit
            end if
         end do
         if (.not. found) then
            nfiles = nfiles + 1
            tmp(nfiles) = fname
         end if
      end do

      if (nfiles > 0) then
         allocate(character(len=maxlen) :: file_paths(nfiles))
         file_paths = tmp(1:nfiles)
      else
         allocate(character(len=1) :: file_paths(0)) 
      end if

      if (allocated(tmp)) deallocate(tmp)
   end subroutine build_file_list

   subroutine scan_and_attach(files, lat0, lon0, y0,mon0,d0,h0,mi0,s0, y1,mon1,d1,h1,mi1,s1, cal_default,        &
                           required_vars_in_file, scans, ok, errmsg)
      character(*), dimension(:), intent(in)  :: files
      real(rk),                   intent(in)  :: lat0, lon0
      integer,                    intent(in)  :: y0,mon0,d0,h0,mi0,s0, y1,mon1,d1,h1,mi1,s1
      character(*),               intent(in)  :: cal_default
      character(*), dimension(:), intent(in)  :: required_vars_in_file   ! Variable names to find in any file
      type(ForcingScan), allocatable, intent(out) :: scans(:)
      logical,                     intent(out) :: ok
      character(*),                intent(out) :: errmsg

      integer :: i
      logical :: ok_inside
      ok = .false.; errmsg = ''      
      allocate(scans(size(files)))
      do i=1, size(files)
         call scan_forcing(files(i), lat0, lon0, y0,mon0,d0,h0,mi0,s0, y1,mon1,d1,h1,mi1,s1,  &
                           cal_default, required_vars_in_file, scans(i), ok_inside, errmsg)
         if (.not. ok_inside) return
      end do
      ok = .true.
   end subroutine scan_and_attach

   subroutine attach_vars_to_scans(vars, files, scans)
      type(VarInfo),              intent(inout) :: vars(:)
      character(*), dimension(:), intent(in)    :: files
      type(ForcingScan), target,  intent(inout) :: scans(:)
      integer :: i, j
      do i=1, size(vars)
         if (associated(vars(i)%scan)) nullify(vars(i)%scan)
         if (len_trim(vars(i)%file_path) == 0) cycle
         do j=1, size(files)
            if (trim(vars(i)%file_path) == trim(files(j))) then
            vars(i)%scan => scans(j)
            vars(i)%calendar  = scans(j)%cal%kind
            vars(i)%first_time = scans(j)%axis%t_first
            vars(i)%last_time  = scans(j)%axis%t_last
            exit
            end if
         end do
         if (.not. associated(vars(i)%scan) .and. vars(i)%input_type == in_file) then
            stop 'Internal: file not found in scans for var '//trim(vars(i)%id)
         end if
      end do
   end subroutine attach_vars_to_scans

   subroutine build_year_indices(var, sim_y_start, sim_y_end, repeat_enabled, repeat_year)
      type(VarInfo), intent(inout) :: var
      integer,       intent(in)    :: sim_y_start, sim_y_end
      logical,       intent(in)    :: repeat_enabled
      integer,       intent(in)    :: repeat_year

      integer :: ny, y, yf, k
      real(rk) :: t0, t1
      ny = sim_y_end - sim_y_start + 1
      if (allocated(var%idx_year)) deallocate(var%idx_year)
      allocate(var%idx_year(2, ny))

      do k=1, ny
         y  = sim_y_start + (k-1)
         yf = merge(repeat_year, y, repeat_enabled)

         ! Build year window in the file’s calendar/units (Jan-01 00:00 to Dec-31 23:59:59 equivalent)
         t0 = seconds_since_datetime_file(var%scan%cal, var%scan%u, yf, 1, 1, 0, 0, 0)
         ! Use next year's start minus a tiny epsilon to be inclusive
         t1 = seconds_since_datetime_file(var%scan%cal, var%scan%u, yf+1, 1, 1, 0, 0, 0)

         ! Clip to the overall available range (and to sim start/end in your driver later if needed)
         if (t0 < var%scan%axis%t_first) t0 = var%scan%axis%t_first
         if (t1 > var%scan%axis%t_last ) t1 = var%scan%axis%t_last

         var%idx_year(1,k) = max(1, index_at_or_after(var%scan%axis%t_s, t0))
         var%idx_year(2,k) = max(var%idx_year(1,k), index_at_or_before(var%scan%axis%t_s, t1))

         ! Add ±1 record padding here if needed -> Check forcing is loaded properly
         ! var%idx_year(1,k) = max(1, var%idx_year(1,k)-1)
         ! var%idx_year(2,k) = min(size(var%scan%axis%t_s), var%idx_year(2,k)+1)
      end do
   end subroutine build_year_indices

   subroutine clear_state(s)
      type(ForcingState), intent(inout) :: s

      ! Reset scalars
      s%cfg%repeat_enabled = .false.
      s%cfg%repeat_year    = -huge(1)
      s%cfg%cfg_calendar   = 0
      s%sim_y_start = 0
      s%sim_y_end   = 0
      s%sim_cal%kind = 0  

      ! Deallocate and re-allocate empty arrays
      if (allocated(s%vars))  deallocate(s%vars)
      if (allocated(s%files)) deallocate(s%files)
      if (allocated(s%scans)) deallocate(s%scans)
      if (allocated(s%scan_idx_of_var)) deallocate(s%scan_idx_of_var)

      allocate(s%vars(0))
      allocate(s%scans(0))
      allocate(character(len=1) :: s%files(0))
      allocate(s%scan_idx_of_var(0))
   end subroutine

   ! Load a single series into a ForcingVarData structure
   subroutine load_var_series(Yv, scan, db, varname, time_name, i0, i1, has_latlon, yi, xi)
      type(ForcingVarData), intent(inout) :: Yv
      type(ForcingScan),    intent(in)    :: scan
      type(NcFile),         intent(in)    :: db
      character(*),         intent(in)    :: varname, time_name
      integer,              intent(in)    :: i0, i1, yi, xi
      logical,              intent(in)    :: has_latlon
      integer :: nt
      ! Auxiliar to calculate time edges
      integer(lk) :: dt_first, dt_last
      integer :: i

      nt = max(0, i1 - i0 + 1)
      if (nt <= 0) stop 'load_var_series: empty time window'

      ! Ensure t_axis has the right size, then fill directly from slice
      if (.not. allocated(Yv%t_axis) .or. size(Yv%t_axis) /= nt) then
         if (allocated(Yv%t_axis)) deallocate(Yv%t_axis)
         allocate(Yv%t_axis(nt))
      end if

      ! Align with simulation time
      Yv%t_axis = int(nint(scan%axis%t_s(i0:i1)), lk) - scan%sim_offset


      ! Ensure data has the right size, then read directly into it
      if (.not. allocated(Yv%data) .or. size(Yv%data) /= nt) then
         if (allocated(Yv%data)) deallocate(Yv%data)
         allocate(Yv%data(nt))
      end if
      call get_ts_slice_from_point( db, trim(varname), trim(time_name), i0, i1, has_latlon, &
                                    trim(scan%lat_name), trim(scan%lon_name), yi, xi, Yv%data )

      ! Cursor metadata
      Yv%is_const = .false.
      Yv%n        = nt

      ! Allocate edges for time
      if (allocated(Yv%t_edge)) deallocate(Yv%t_edge)
      allocate(Yv%t_edge(Yv%n + 1))

      if (Yv%n >= 2) then
         dt_first = max(1_lk, Yv%t_axis(2) - Yv%t_axis(1))
         dt_last  = max(1_lk, Yv%t_axis(Yv%n) - Yv%t_axis(Yv%n-1))

         ! Extrapolated outer edges
         Yv%t_edge(1)       = Yv%t_axis(1) - dt_first/2_lk
         Yv%t_edge(Yv%n+1)  = Yv%t_axis(Yv%n) + dt_last/2_lk

         ! Midpoints between consecutive stamps
         do i = 1, Yv%n-1
            Yv%t_edge(i+1) = (Yv%t_axis(i) + Yv%t_axis(i+1)) / 2_lk  ! integer midpoint
         end do
      else
         ! n == 1: make a very wide interval so it never advances
         Yv%t_edge(1)      = Yv%t_axis(1) - huge(1_lk)/4_lk
         Yv%t_edge(2)      = Yv%t_axis(1) + huge(1_lk)/4_lk
      end if
      Yv%idx      = 1
      Yv%t_next   = merge(Yv%t_axis(2), huge(1_lk), nt >= 2)
   end subroutine load_var_series



   pure integer function find_var_index(vars, id) result(idx)
      type(VarInfo), intent(in) :: vars(:)
      character(*),  intent(in) :: id
      integer :: k
      idx = 0
      do k=1, size(vars)
         if (trim(vars(k)%id) == trim(id)) then; idx = k; return; end if
      end do
   end function find_var_index

   pure integer function file_index_of_var(v, scans) result(j)
      type(VarInfo),     intent(in) :: v
      type(ForcingScan), intent(in) :: scans(:)
      integer :: k
      j = 0
      do k=1, size(scans)
         if (associated(v%scan) .and. (v%scan%time_name == scans(k)%time_name)) then
            j = k; return
         end if
      end do
      ! Fallback if time_name isn't unique; compare pointer address by file fields if you store them,
      ! or pass 'j' from the caller. For now assume time_name unique per file in your dataset.
   end function file_index_of_var

end module load_forcing
