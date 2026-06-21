module event_types
    use precision_types, only: rk

    implicit none
    private

    public :: Event, EventData
    public :: EVENT_NONE, EVENT_TRACER_PULSE, EVENT_TRAWLING

    character(len=*), parameter :: EVENT_NONE         = 'none'
    character(len=*), parameter :: EVENT_TRACER_PULSE = 'tracer_pulse'
    character(len=*), parameter :: EVENT_TRAWLING     = 'trawling'

    ! Polymorphic type to store any event details 
    ! This will be extended by the particular event module
    type, abstract :: EventData
    end type EventData
    
    type :: Event
        character(len=:), allocatable :: name
        character(len=:), allocatable :: event_type
        character(len=:), allocatable :: date_string
        ! Time
        real(rk) :: t_start  = 0._rk
        real(rk) :: t_end    = 0._rk
        real(rk) :: duration = 0._rk
        integer  :: handler_idx = 0
        integer  :: definition_idx = 0
        integer  :: occurrence_idx = 0
        logical  :: is_instantaneous = .false.
        logical  :: changes_tendency = .false.
        logical  :: was_triggered    = .false.
        ! Specific details for the event
        class(EventData), allocatable :: details
    contains
        procedure :: clear          => event_clear
        procedure :: is_active      => event_is_active
        procedure :: should_trigger => event_should_trigger
    end type Event

contains

    logical function event_is_active(self, model_time)
        class(Event), intent(in) :: self
        real(rk),     intent(in) :: model_time

        event_is_active = .false.

        if (self%is_instantaneous) return

        event_is_active = (model_time >= self%t_start) .and. &
                          (model_time <  self%t_end)
    end function event_is_active

    logical function event_should_trigger(self, model_time)
        class(Event), intent(in) :: self
        real(rk),     intent(in) :: model_time

        real(rk) :: eps_t

        event_should_trigger = .false.

        if (.not. self%is_instantaneous) return
        if (self%was_triggered) return

        eps_t = 1.0e-9_rk * max(1.0_rk, abs(model_time), abs(self%t_start))

        event_should_trigger = abs(model_time - self%t_start) <= eps_t
    end function event_should_trigger

    subroutine event_clear(self)
        class(Event), intent(inout) :: self

        if (allocated(self%name))        deallocate(self%name)
        if (allocated(self%event_type))  deallocate(self%event_type)
        if (allocated(self%date_string)) deallocate(self%date_string)
        if (allocated(self%details))     deallocate(self%details)

        self%t_start     = 0._rk
        self%t_end       = 0._rk
        self%duration    = 0._rk
        self%handler_idx = 0
        self%definition_idx    = 0
        self%occurrence_idx    = 0
        self%is_instantaneous  = .false.
        self%was_triggered     = .false.
        self%changes_tendency  = .false.
    end subroutine event_clear

end module event_types