module data_types
   use precision_types, only: rk, lk
   use time_types,      only: CFCalendar, CFUnits
   implicit none

   private

   public :: DATA_INPUT_FILE, DATA_INPUT_CONSTANT, DATA_INPUT_COMPUTE, DATA_INPUT_OFF
   public :: DATA_FORMAT_UNKNOWN, DATA_FORMAT_NETCDF, DATA_FORMAT_CSV
   public :: DATA_TIME_ABSOLUTE, DATA_TIME_REPEAT_YEAR
   public :: TimeWindow
   public :: DataLoaderCfg
   public :: DataSpec
   public :: DataFileInfo
   public :: DataVarSeries
   public :: InputData

   integer, parameter :: DATA_INPUT_FILE     = 1
   integer, parameter :: DATA_INPUT_CONSTANT = 2
   integer, parameter :: DATA_INPUT_COMPUTE  = 3
   integer, parameter :: DATA_INPUT_OFF      = 4

   integer, parameter :: DATA_FORMAT_UNKNOWN = 0
   integer, parameter :: DATA_FORMAT_NETCDF  = 1
   integer, parameter :: DATA_FORMAT_CSV     = 2

   integer, parameter :: DATA_TIME_ABSOLUTE    = 1
   integer, parameter :: DATA_TIME_REPEAT_YEAR = 2

   type :: TimeWindow
      real(rk) :: t_start = 0.0_rk
      real(rk) :: t_end   = 0.0_rk
   end type TimeWindow

   type :: DataLoaderCfg
      logical :: load_yearly    = .true.
      logical :: repeat_enabled = .false.
      integer :: repeat_year    = -huge(1)
      integer :: cfg_calendar   = 0
      integer :: time_mode      = DATA_TIME_ABSOLUTE
   end type DataLoaderCfg

   type :: DataSpec
      character(:), allocatable :: name        ! internal/model name
      character(:), allocatable :: source_var  ! variable name in file
      character(:), allocatable :: path        ! source file path
      character(:), allocatable :: time_var    ! e.g. "time" or "timestamp"
      character(:), allocatable :: units

      integer :: input_type = DATA_INPUT_OFF
      integer :: format     = DATA_FORMAT_UNKNOWN

      logical  :: has_location = .false.
      real(rk) :: lon = 0.0_rk
      real(rk) :: lat = 0.0_rk

      real(rk) :: const_value = 0.0_rk

      integer :: file_index = 0
      integer :: calendar   = 0

      real(rk) :: first_time = 0.0_rk
      real(rk) :: last_time  = 0.0_rk

      integer, allocatable :: idx_window(:,:)
   end type DataSpec

   type :: DataFileInfo
      character(:), allocatable :: path
      character(:), allocatable :: extension
      integer :: format = DATA_FORMAT_UNKNOWN
   end type DataFileInfo

   type :: DataVarSeries
      character(:), allocatable :: name
      character(:), allocatable :: units

      logical  :: is_const = .false.
      real(rk) :: const_value = 0.0_rk

      integer(lk), allocatable :: t_axis(:)
      integer(lk), allocatable :: t_edge(:)
      real(rk),    allocatable :: values(:)

      type(CFCalendar) :: cal
      type(CFUnits)    :: u
      integer(lk)      :: sim_offset = 0_lk
      integer          :: time_mode = DATA_TIME_ABSOLUTE
      integer          :: repeat_year = -1

      integer     :: idx = 1
      integer     :: n   = 0
      integer(lk) :: t_next = huge(1_lk)
   end type DataVarSeries

   type :: InputData
      type(DataVarSeries), allocatable :: vars(:)
   end type InputData

end module data_types