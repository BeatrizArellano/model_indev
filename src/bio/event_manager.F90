module event_manager
   use fabm,             only: type_fabm_model
   use grids,            only: VerticalGrid
   use precision_types,  only: rk
   use read_config_yaml, only: ConfigParams
   use time_types,       only: DateTime, CFCalendar
   use tracer_pulse,     only: TracerPulse, tracer_pulse_init, tracer_pulse_apply, tracer_pulse_clear


   implicit none
   private

   public :: EventManager

   type :: EventManager
      logical :: is_initialized =   .false.
      logical :: any_events     =   .false.
      logical :: require_tendency = .false.
      !--- Tracer Pulse ---
      type(TracerPulse) :: TP
      logical :: tracer_pulse_enabled = .false.
      character(len=:), allocatable :: tracer_pulse_cfgfile
      type(ConfigParams) :: tracer_pulse_cfg
   contains
      procedure :: init  => em_init
      procedure :: clear => em_clear
      procedure :: apply => em_apply
   end type EventManager

contains

   subroutine em_init(self, main_cfg, grid, FabmMod, sim_startdate, sim_enddate, cal)
      class(EventManager),    intent(inout) :: self
      type(ConfigParams),     intent(in)    :: main_cfg
      type(VerticalGrid),     intent(in)    :: grid
      class(type_fabm_model), pointer, intent(in)    :: FabmMod
      type(DateTime),         intent(in)    :: sim_startdate, sim_enddate
      type(CFCalendar),       intent(in)    :: cal

      character(len=256) :: key

      call self%clear()

      ! -------------------------------
      ! tracer_pulse activation & config
      ! -------------------------------
      key = 'biogeochemistry.process_modules.tracer_pulse.enabled'
      if (main_cfg%is_set(key)) then
         self%tracer_pulse_enabled = main_cfg%get_param_logical(key, default=.false., strict=.false.)
      else
         self%tracer_pulse_enabled = .false.
      end if

      if (self%tracer_pulse_enabled) then
         key = 'biogeochemistry.process_modules.tracer_pulse.config_file'

         if (.not. main_cfg%is_set(key)) then
            call fatal('event_manager:init', 'tracer_pulse enabled but "'//key//'" missing/null.')
         end if

         self%require_tendency = .true. 

         self%tracer_pulse_cfgfile = main_cfg%get_param_str(key, required=.true., trim_value=.true.)

         call self%tracer_pulse_cfg%init()
         call self%tracer_pulse_cfg%load_yaml_content(self%tracer_pulse_cfgfile)

         ! Initialise tracer pulse
         call tracer_pulse_init(self%tracer_pulse_cfg, grid, FabmMod, sim_startdate, sim_enddate, cal, self%TP)
      end if

      self%any_events = self%tracer_pulse_enabled

      self%is_initialized = .true.
   end subroutine em_init


   subroutine em_apply(self, model_time, dt, tendency)
      class(EventManager), intent(inout) :: self
      real(rk),            intent(in)    :: model_time, dt
      real(rk),            intent(inout) :: tendency(:,:)   ! Tendency array for all tracers

      if (.not. self%is_initialized) then
         call fatal('event_manager:apply', 'EventManager not initialized.')
      end if

      ! For now: no-op. Next step: call modules and add tendencies.
      if (self%tracer_pulse_enabled) then
         call tracer_pulse_apply(model_time, dt, self%TP, tendency)
      end if
   end subroutine em_apply


   subroutine em_clear(self)
      class(EventManager), intent(inout) :: self
      self%tracer_pulse_enabled = .false.
      if (allocated(self%tracer_pulse_cfgfile)) deallocate(self%tracer_pulse_cfgfile)
      call self%tracer_pulse_cfg%clear()
      call tracer_pulse_clear(self%TP)
      self%is_initialized = .false.
      self%any_events       = .false.
      self%require_tendency = .false.
   end subroutine em_clear


   ! ---- Local error routine ----
   subroutine fatal(where, msg)
      character(len=*), intent(in) :: where, msg
      write(*,'(A,1X,A)') '[FATAL:'//trim(where)//']', trim(msg)
      stop 1
   end subroutine fatal

end module event_manager
