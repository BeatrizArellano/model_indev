!===============================================================================
!  Module: read_config_yaml
!  Reads and loads parameters from a YAML config file using FABM's yaml parsing libraries.
!  - Supports full dotted parameters ("time.dt") and unique name-only parameters ("dt")
!  - Numeric: returns real(rk); parses ints/reals/boolean-like strings
!  - Logical: true/false/yes/no/on/off/0/1 (case-insensitive)
!  Quick start:
!    use read_config_yaml, only: ConfigParams
!    type(ConfigParams) :: main_cfg, physics_cfg
!    call main_cfg%init();    call main_cfg%load_yaml_content('main.yaml')
!    call physics_cfg%init(); call physics_cfg%load_yaml_content('physics.yaml')
!
!  Scalars with validation
!    real(rk)    :: dt    = main_cfg%get_param_num('time.dt', required=.true., positive=.true.)
!    character(:), allocatable :: ofile
!    ofile = main_cfg%get_param_str('output.file', default='out.nc', trim_value=.true.)
!
!  Booleans (strict or relaxed)
!    logical :: from_file
!    from_file = physics_cfg%get_param_logical('tides.read_from_file', required=.true.)         ! strict: true/false
!    from_file = physics_cfg%get_param_logical('tides.read_from_file', required=.true., strict=.false.) ! relaxed: yes/no,on/off,1/0
!
!  Null-handling & overrides
!    if (.not. physics_cfg%is_null('forcing.surf_air_temp.constant')) then
!       real(rk) :: c = physics_cfg%get_param_num('forcing.surf_air_temp.constant', finite=.true.)
!       ! use constant override
!    else
!       character(:), allocatable :: vfile
!       if (physics_cfg%is_set('forcing.surf_air_temp.filename')) then
!          vfile = physics_cfg%get_param_str('forcing.surf_air_temp.filename', trim_value=.true.)
!       else
!          vfile = physics_cfg%get_param_str('forcing.filename', required=.true.)
!       end if
!       ! read from vfile/name
!    end if
!  ConfigParams Object:
!    * Each ConfigParams instance is independent (use one per YAML file).
!    * Reuse: call %clear() to drop previous configuration trees,
!             then %load_yaml_content(...) again.
!    * End clearing resources using call %clear().
!
!  YAML conventions:
!    * Null / “unset”:  ~, null/NULL/Null, none/None, no, "", or empty scalar (key:)
!      - Missing and null are treated the same by getters and predicates.
!    * Booleans:
!        - Preferred: true/false (case-insensitive).
!        - Relaxed mode (strict=.false.): yes/no, on/off, 1/0 (case-insensitive).
!    * Numbers: integers, realsare accepted (e.g., 3, 3.0).
!
!  API overview:
!    Predicates:
!      logical :: has(path)              ! present in any level (incl. null)
!      logical :: is_absent(path)        ! not present 
!      logical :: is_null(path[, include_empty_scalar, coerce_none])
!      logical :: is_set(path[, include_empty_scalar, coerce_none]) ! present and not null
!      logical :: is_numeric(path)
!      logical :: is_boolean(path[, relaxed])
!      logical :: is_string(path[, include_empty_scalar])
!
!    Getters:
!      real(rk)    function get_param_num(path,    &
!                          default, found, required, min, max,        &
!                          positive, nonnegative, finite)
!      integer     function get_param_int(path,    &
!                          default, found, required, min, max)        
!      logical     function get_param_logical(path,                    &
!                          default, found, required, strict)           ! strict=.true. by default
!      character(:), allocatable function get_param_str(path,          &
!                          default, found, required, choices,          &
!                          trim_value, match_case, empty_ok, allow_numeric)
!
!    Utilities:
!      logical :: has_param(path)      ! synonym of has(path)
!      character(len=:), allocatable :: list(:) = list_params(prefix)  ! includes null keys
!   Behaviour:
!    * Required vs default:
!        - required=.true. -> missing or null raises a clear error.
!        - default=...     -> returned if missing or null; found=.false.
!    * Validation (numbers): positive / nonnegative / min / max / finite (NaN check).
!    * Strings: trimming by default; optional enum check via choices=(/ 'Euler','Other' /)
!    * Name-only lookup: if a short parameter is unique across levels, it resolves without a
!      full dotted path; ambiguous short names return an error.
!
!=====================================================================================================
module read_config_yaml
   use, intrinsic :: ieee_arithmetic
   use precision_types, only: rk   
   use yaml_types
   use yaml, only: yaml_parse => parse, yaml_error_length => error_length
   use str_utils, only: to_lower

   implicit none
   private
 
   public :: ConfigParams

   integer, parameter :: PARAMLEN = 256
   integer, parameter :: STRLEN = 512

   type cfg_real_t
      character(len=PARAMLEN) :: key = ''
      real(rk)              :: val = 0.0_rk
   end type
   type cfg_int_t
      character(len=PARAMLEN) :: key = ''
      integer               :: val = 0
   end type
   type cfg_logical_t
      character(len=PARAMLEN) :: key = ''
      logical               :: val = .false.
   end type
   type cfg_string_t
      character(len=PARAMLEN) :: key = ''
      character(len=STRLEN) :: val = ''
   end type

   type cfg_null_t
      character(len=PARAMLEN) :: key = ''
   end type


   type :: ConfigParams
      private
      type(cfg_real_t),    allocatable :: store_reals(:)
      type(cfg_int_t),     allocatable :: store_ints(:)
      type(cfg_logical_t), allocatable :: store_logs(:)
      type(cfg_string_t),  allocatable :: store_strs(:)
      type(cfg_null_t),    allocatable :: store_nulls(:)
   contains
      procedure :: init                => cfg_init
      procedure :: clear               => cfg_clear
      procedure :: load_yaml_content   => cfg_load_yaml_content
      procedure :: has_param           => cfg_has_param
      procedure :: list_params         => cfg_list_params
      procedure :: is_absent           => cfg_is_absent
      procedure :: is_null             => cfg_is_null
      procedure :: is_disabled         => cfg_is_disabled
      procedure :: is_set              => cfg_is_set
      procedure :: is_numeric          => cfg_is_numeric
      procedure :: is_boolean          => cfg_is_boolean
      procedure :: is_string           => cfg_is_string
      procedure :: get_param_str       => cfg_get_param_str
      procedure :: get_param_num       => cfg_get_param_num
      procedure :: get_param_int       => cfg_get_param_int
      procedure :: get_param_logical   => cfg_get_param_logical
      ! internals
      procedure, private :: walk_dictionary
      procedure, private :: store_scalar
      procedure, private :: resolve_key
      procedure, private :: exists_fullkey    ! Once the full key to the parameter has been found
      procedure, private :: lookup_string
      procedure, private :: append_real
      procedure, private :: append_int
      procedure, private :: append_log
      procedure, private :: append_str
      procedure, private :: find_real_key
      procedure, private :: find_int_key
      procedure, private :: find_log_key
      procedure, private :: find_str_key
      procedure, private :: append_null
      procedure, private :: find_null_key
   end type ConfigParams

contains
   subroutine cfg_init(self)
      class(ConfigParams), intent(inout) :: self
      call self%clear()
   end subroutine cfg_init

   subroutine cfg_clear(self)
      class(ConfigParams), intent(inout) :: self
      if (allocated(self%store_reals)) deallocate(self%store_reals)
      if (allocated(self%store_ints))  deallocate(self%store_ints)
      if (allocated(self%store_logs))  deallocate(self%store_logs)
      if (allocated(self%store_strs))  deallocate(self%store_strs)
      if (allocated(self%store_nulls)) deallocate(self%store_nulls) 
   end subroutine cfg_clear

   subroutine cfg_load_yaml_content(self, path, prefix)
      class(ConfigParams), intent(inout) :: self
      character(len=*),  intent(in)     :: path
      character(len=*),  intent(in), optional :: prefix
      class(type_node),       pointer :: root
      class(type_dictionary), pointer :: dict
      character(len=yaml_error_length) :: err
      integer :: unit
      character(len=PARAMLEN) :: pre

      pre = ''; if (present(prefix)) pre = trim(prefix)
      unit = 99

      root => yaml_parse(trim(path), unit, err)
      if (err /= '') call stop_missingpar('load_yaml_content', trim(err))
      if (.not.associated(root)) call stop_missingpar('load_yaml_content', 'Empty YAML: '//trim(path))

      select type (root)
      class is (type_dictionary)
         dict => root
      class default
         call stop_missingpar('load_yaml_content', trim(path)//' must have a dictionary at root.')
      end select

      call self%walk_dictionary(dict, trim(pre))
   end subroutine cfg_load_yaml_content

   recursive subroutine walk_dictionary(self, mapping, prefix)
      class(ConfigParams),      intent(inout) :: self
      class(type_dictionary),   intent(in)    :: mapping
      character(len=*),         intent(in)    :: prefix
      type(type_key_value_pair), pointer :: pair
      class(type_node), pointer           :: node
      character(len=PARAMLEN) :: joined

      pair => mapping%first
      do while (associated(pair))
         node => pair%value
         joined = join_key(prefix, pair%key)
         select type (node)
         class is (type_dictionary)
            call self%walk_dictionary(node, trim(joined))
         class is (type_scalar)
            call self%store_scalar(trim(joined), node)
         class is (type_null)
            call self%append_null(trim(joined))
         class default
            ! sequences/others: ignore for now, can be added later
         end select
         pair => pair%next
      end do
   end subroutine walk_dictionary

   subroutine store_scalar(self, flat_key, s)
      class(ConfigParams), intent(inout) :: self
      character(len=*),    intent(in)    :: flat_key
      class(type_scalar),  intent(in)    :: s
      logical :: ok
      real(rk) :: rval
      integer  :: ios, ival
      character(len=STRLEN) :: str, low

      ! Treat empty / "null" / "~" / "none"/'nil' as YAML null ---
      str = trim(adjustl(s%string))
      low = to_lower(str)
      if (len_trim(str) == 0 .or. low=='null' .or. low=='~' .or. low=='none' .or. low=='nil') then
         call self%append_null(flat_key)
         return
      end if

      ! Try real (covers ints in many builds too)
      rval = s%to_real(default=0.0_rk, success=ok)
      if (ok) then
         call self%append_real(flat_key, rval); return
      end if

      ! Try parse integer from string
      str = adjustl(trim(s%string))
      read(str, *, iostat=ios) ival
      if (ios == 0) then
         call self%append_int(flat_key, ival); return
      end if

      ! Try boolean-like strings
      low = to_lower(str)
      if (is_true_string(low))  then; call self%append_log(flat_key, .true. ); return; end if
      if (is_false_string(low)) then; call self%append_log(flat_key, .false.); return; end if

      ! Fallback: store as string
      call self%append_str(flat_key, str)
   end subroutine store_scalar

   !==================== Getters ====================

   function cfg_get_param_str(self, key, default, found, required, choices, trim_value, match_case, empty_ok, allow_numeric) result(s)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      character(len=*),    intent(in),  optional :: default
      logical,             intent(out), optional :: found
      logical,             intent(in),  optional :: required
      character(len=*),    intent(in),  optional :: choices(:)
      logical,             intent(in),  optional :: trim_value     ! default: .true.
      logical,             intent(in),  optional :: match_case     ! default: .false.
      logical,             intent(in),  optional :: empty_ok       ! default: .false.
      logical,             intent(in),  optional :: allow_numeric  ! default: .true.

      character(len=STRLEN) :: s
      character(len=PARAMLEN) :: full
      logical :: req, do_trim, case_sens, ok_empty, coerce, match
      integer :: i, j
      character(len=STRLEN) :: tmp, s_lc, c_lc

      full      = self%resolve_key(key)
      req       = .false.; if (present(required))    req      = required
      do_trim   = .true.;  if (present(trim_value))  do_trim  = trim_value
      case_sens = .false.; if (present(match_case)) case_sens  = match_case
      ok_empty  = .false.; if (present(empty_ok))    ok_empty = empty_ok
      coerce    = .true.;  if (present(allow_numeric)) coerce  = allow_numeric

      ! Check if it's set (the parameter is present and is not null) (treat YAML null/~ as "not set")
      if ( self%is_absent(key) .or. self%is_null(key) ) then
         if (present(found)) found = .false.
         if (present(default)) then
            s = default
            if (do_trim) s = trim(adjustl(s))
            if (present(found)) found = .false.
            return
         else if (req) then
            call stop_missingpar('get_param_str','Required parameter "'//trim(key)//'" is missing or null.')
         else
            call stop_missingpar('get_param_str','Parameter "'//trim(key)//'" is missing or null and no default provided.')
         end if
      end if

      ! Read string
      i = self%find_str_key(full)
      if (i /= 0) then
         s = self%store_strs(i)%val
      else if (coerce) then
         ! Convert numeric/logical  to string if allowed
         if (self%find_real_key(full) /= 0) then
            write(s,'(G0.12)') self%store_reals(self%find_real_key(full))%val
         else if (self%find_int_key(full) /= 0) then
            write(s,'(I0)') self%store_ints(self%find_int_key(full))%val
         else if (self%find_log_key(full) /= 0) then
            if (self%store_logs(self%find_log_key(full))%val) then
               s = 'true'
            else
               s = 'false'
            end if
         else
            ! Try lookup_string (e.g., numeric as text)
            tmp = adjustl(trim(self%lookup_string(full)))
            if (len_trim(tmp) == 0) then
               ! Treat as Empty string
               s = ''
            else
               s = tmp
            end if
         end if
      else
         call stop_missingpar('get_param_str','Parameter "'//trim(key)//'" is not a string.')
      end if

      ! Trim
      if (do_trim) s = trim(adjustl(s))

      ! Empty-string policy
      if (len_trim(s) == 0 .and. .not. ok_empty) then
         if (present(default)) then
            s = default
            if (do_trim) s = trim(adjustl(s))
            if (present(found)) found = .false.
            return
         else
            call stop_missingpar('get_param_str','"'//trim(key)//'" must be a non-empty string.')
         end if
      end if

      ! Check for choices if present
      if (present(choices)) then
         match = .false.
         if (.not. case_sens) then            
            s_lc = to_lower(s)
            do j = 1, size(choices)
               c_lc = to_lower( trim(choices(j)) )
               if (s_lc == c_lc) then
                  match = .true.; exit
               end if
            end do
         else            
            do j = 1, size(choices)
               if (s == trim(choices(j))) then
                  match = .true.; exit
               end if
            end do
         end if
         if (.not. match) then
            call stop_missingpar('get_param_str','"'//trim(key)//'" must be one of the allowed choices.')
         end if
      end if

      if (present(found)) found = .true.
   end function cfg_get_param_str


   function cfg_get_param_logical(self, key, default, found, required, strict) result(val)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      logical, optional,   intent(in) :: default, required, strict
      logical, optional,   intent(out):: found
      logical :: val, req, sstrict
      character(len=PARAMLEN) :: full
      integer :: i
      character(len=STRLEN) :: s

      full = self%resolve_key(key)
      req  = .false.; if (present(required)) req = required
      sstrict = .false.; if (present(strict)) sstrict = strict

      ! Detecting if missing or null
      if ( self%is_absent(key) .or. self%is_null(key) ) then
         if (present(found)) found = .false.
         if (present(default)) then
            val = default; return
         else if (req) then
            call stop_missingpar('get_param_logical','"'//trim(key)//'" is missing or null.')
         else
            call stop_missingpar('get_param_logical','"'//trim(key)//'" missing/null and no default.')
         end if
      end if

      ! Check whether it is boolean
      if (.not. self%is_boolean(key, relaxed=.not. sstrict)) then
         call stop_missingpar('get_param_logical','Parameter "'//trim(key)//'" is not boolean.')
      end if

      ! Extract
      i = self%find_log_key(full)
      if (i /= 0) then
         val = self%store_logs(i)%val
      else
         s = trim(adjustl(self%lookup_string(full)))
         s = to_lower(s)
         if (s == 'true') then
            val = .true.
         else if (s == 'false') then
            val = .false.
         else if (.not. sstrict) then
            ! Potential alternatives (synonims)
            if (s=='yes' .or. s=='on'  .or. s=='1') then
               val = .true.
            else if (s=='no'  .or. s=='off' .or. s=='0') then
               val = .false.
            else
               call stop_missingpar('get_param_log','Unrecognized boolean literal "'//trim(s)//'".')
            end if
         else
            call stop_missingpar('get_param_log','Unrecognized boolean literal "'//trim(s)//'".')
         end if
      end if

      if (present(found)) found = .true.
   end function cfg_get_param_logical


   function cfg_get_param_num(self, key, default, found, required, min, max, positive, nonnegative, finite) result(val)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      real(rk), optional, intent(in) :: default
      logical,  optional, intent(out):: found
      logical,  optional, intent(in) :: required, positive, nonnegative, finite
      real(rk), optional, intent(in) :: min, max

      logical :: req
      real(rk) :: val
      character(len=PARAMLEN) :: full
      character(len=STRLEN) :: s
      integer :: i, ios
      real(rk) :: rtmp

      full = self%resolve_key(key)
      req = .false.; if (present(required)) req = required     

      if (.not. self%is_set(key)) then
         if (present(found)) found = .false.
         if (present(default)) then
            val = default
            return
         else if (req) then
            call stop_missingpar('get_param_num','Required parameter "'//trim(key)//'" is missing or null.')
         else
            ! no default, not required: choose a policy; safest is to error
            call stop_missingpar('get_param_num','Parameter "'//trim(key)//'" is missing or null and no default provided.')
         end if
      end if

      ! Stop if not nuemric
      if (.not. self%is_numeric(key)) then
         call stop_missingpar('get_param_num','Parameter "'//trim(key)//'" is not numeric.')
      end if

      ! 2) Extract value (prefer stored numeric; else parse numeric string)
      i = self%find_real_key(full)
      if (i /= 0) then
         val = self%store_reals(i)%val
      else
         i = self%find_int_key(full)
         if (i /= 0) then
            val = real(self%store_ints(i)%val, rk)
         else
            s = adjustl(trim(self%lookup_string(full)))
            read(s, *, iostat=ios) val
            if (ios /= 0) then
               ! (fallback for compilers that don’t read reals directly first)
               read(s, *, iostat=ios) rtmp
               if (ios == 0) val = rtmp
            end if
         end if
      end if

      ! Validate 
      if (present(positive)     .and.  positive     .and. .not.(val > 0.0_rk))  &
         call stop_missingpar('get_param_num','"'//trim(key)//'" must be > 0.')
      if (present(nonnegative)  .and.  nonnegative  .and.      (val < 0.0_rk))  &
         call stop_missingpar('get_param_num','"'//trim(key)//'" must be >= 0.')
      if (present(min)          .and.               (val < min))                &
         call stop_missingpar('get_param_num', limit_msg_num(key, '>=', min, val))
      if (present(max)          .and.               (val > max))                &
         call stop_missingpar('get_param_num', limit_msg_num(key, '<=', max, val))
      if (present(finite) .and. finite .and. ieee_is_nan(val))    &
         call stop_missingpar('get_param_num','"'//trim(key)//'" must be finite.')

      if (present(found)) found = .true.
   end function cfg_get_param_num


   integer function cfg_get_param_int(self, key, default, found, required, min, max, positive) result(ival)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      integer,  optional,  intent(in) :: default, min, max
      logical,  optional,  intent(out):: found
      logical,  optional,  intent(in) :: required, positive
      real(rk) :: r
      logical :: has
      integer :: dmin, dmax

      ! get as real with suitable bounds
      if (present(default)) then
         r = self%get_param_num(key, default=real(default, rk), found=has, finite=.true.)
      else
         r = self%get_param_num(key, required=merge(.true., .false., present(required)), found=has, finite=.true.)
      end if

      ! integerness check
      if (abs(r - nint(r)) > 1.0e-6_rk) then
         call stop_missingpar('get_param_int','"'//trim(key)//'" must be an integer.')
      end if
      ival = nint(r)

      ! integer range check (after cast)
      if (present(min)) then
         if (ival < min) call stop_missingpar('get_param_int','"'//trim(key)//'" must be >= min.')
      end if
      if (present(max)) then
         if (ival > max) call stop_missingpar('get_param_int','"'//trim(key)//'" must be <= max.')
      end if
      if (present(positive)     .and.  positive     .and. .not.(ival > 0))  then
         call stop_missingpar('get_param_int','"'//trim(key)//'" must be > 0.')
      end if

      if (present(found)) found = has
   end function cfg_get_param_int


   !==================== Queries ====================

   function cfg_has_param(self, key) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      logical :: tf
      character(len=PARAMLEN) :: full
      full = self%resolve_key(key)
      if (len_trim(full) == 0) then
         tf = .false.
      else
         tf = (self%find_real_key(full) /= 0) .or. (self%find_int_key(full) /= 0) .or. &
            (self%find_log_key(full)  /= 0) .or. (self%find_str_key(full)  /= 0) .or. &
            (self%find_null_key(full) /= 0)
      end if
   end function cfg_has_param


   function cfg_list_params(self, prefix) result(list)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: prefix
      character(len=PARAMLEN), allocatable :: list(:)
      integer :: n, i
      character(len=PARAMLEN) :: pre

      pre = trim(prefix)
      n = 0
      if (allocated(self%store_reals)) then
         do i=1,size(self%store_reals); if (startswith(self%store_reals(i)%key, pre)) n=n+1; end do
      end if
      if (allocated(self%store_ints)) then
         do i=1,size(self%store_ints); if (startswith(self%store_ints(i)%key, pre)) n=n+1; end do
      end if
      if (allocated(self%store_logs)) then
         do i=1,size(self%store_logs); if (startswith(self%store_logs(i)%key, pre)) n=n+1; end do
      end if
      if (allocated(self%store_strs)) then
         do i=1,size(self%store_strs); if (startswith(self%store_strs(i)%key, pre)) n=n+1; end do
      end if
      if (allocated(self%store_nulls)) then
         do i=1,size(self%store_nulls); if (startswith(self%store_nulls(i)%key, pre)) n=n+1; end do
      end if

      allocate(list(n))
      n = 0
      if (allocated(self%store_reals)) then
         do i=1,size(self%store_reals)
            if (startswith(self%store_reals(i)%key, pre)) then; n=n+1; list(n)=self%store_reals(i)%key; end if
         end do
      end if
      if (allocated(self%store_ints)) then
         do i=1,size(self%store_ints)
            if (startswith(self%store_ints(i)%key, pre)) then; n=n+1; list(n)=self%store_ints(i)%key; end if
         end do
      end if
      if (allocated(self%store_logs)) then
         do i=1,size(self%store_logs)
            if (startswith(self%store_logs(i)%key, pre)) then; n=n+1; list(n)=self%store_logs(i)%key; end if
         end do
      end if
      if (allocated(self%store_strs)) then
         do i=1,size(self%store_strs)
            if (startswith(self%store_strs(i)%key, pre)) then; n=n+1; list(n)=self%store_strs(i)%key; end if
         end do
      end if
      if (allocated(self%store_nulls)) then
         do i=1,size(self%store_nulls)
            if (startswith(self%store_nulls(i)%key, pre)) then
               n=n+1; list(n)=self%store_nulls(i)%key
            end if
         end do
      end if
   end function cfg_list_params
   
   logical function cfg_is_absent(self, key) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      character(len=PARAMLEN) :: full
      full = self%resolve_key(key)
      tf = (len_trim(full) == 0)
   end function cfg_is_absent

   logical function cfg_is_null(self, key, include_empty_scalar, coerce_none) result(isn)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      logical, intent(in), optional :: include_empty_scalar   ! default .true.
      logical, intent(in), optional :: coerce_none            ! default .false.
      character(len=PARAMLEN) :: full
      logical :: incl, coer
      character(len=STRLEN) :: s
      integer :: i

      incl = .true.; if (present(include_empty_scalar)) incl = include_empty_scalar
      coer = .true.; if (present(coerce_none)) coer = coerce_none

      isn = .false.
      full = self%resolve_key(key)
      if (len_trim(full) == 0) return

      ! 1) explicit YAML null
      if (self%find_null_key(full) /= 0) then
         isn = .true.; return
      end if

      ! 2) treat empty scalar or "None" as null if requested
      i = self%find_str_key(full)
      if (i /= 0) then
         s = trim(self%store_strs(i)%val)
         if (incl .and. len_trim(s) == 0) then
            isn = .true.; return
         end if
         if (coer .and. (s == 'None' .or. s == 'none' .or. s == 'NONE')) then
            isn = .true.; return
         end if
      end if
   end function cfg_is_null

   logical function cfg_is_set(self, key, include_empty_scalar, coerce_none) result(ok)
         class(ConfigParams), intent(in) :: self
         character(len=*),   intent(in) :: key
         logical, intent(in), optional :: include_empty_scalar, coerce_none
         ok = self%has_param(key) .and. .not. self%is_null(key, include_empty_scalar, coerce_none)
      end function cfg_is_set

      logical function cfg_is_numeric(self, key) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      character(len=PARAMLEN) :: full
      character(len=STRLEN)   :: s
      integer :: i, ios, itmp
      real(rk) :: rtmp

      tf = .false.
      full = self%resolve_key(key)
      if (len_trim(full) == 0) return
      if (self%is_null(key)) return

      ! Direct numeric stores first
      if (self%find_real_key(full) /= 0) then; tf = .true.; return; end if
      if (self%find_int_key(full)  /= 0) then; tf = .true.; return; end if

      ! Booleans are NOT numeric:
      if (self%find_log_key(full)  /= 0) return

      ! Try parse from string store
      i = self%find_str_key(full)
      if (i == 0) return
      s = adjustl(trim(self%store_strs(i)%val))

      read(s, *, iostat=ios) itmp; if (ios == 0) then; tf = .true.; return; end if
      read(s, *, iostat=ios) rtmp; if (ios == 0) then; tf = .true.; return; end if
   end function cfg_is_numeric

   logical function cfg_is_boolean(self, key, relaxed) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      logical, optional,   intent(in) :: relaxed
      logical :: rel
      character(len=PARAMLEN) :: full
      tf = .false.
      rel = .false.; if (present(relaxed)) rel = relaxed

      full = self%resolve_key(key)
      if (len_trim(full) == 0) return
      if (self%is_null(key)) return

      if (self%find_log_key(full) /= 0) then
         tf = .true.; return
      end if

      block
         integer :: i
         character(len=STRLEN) :: s
         i = self%find_str_key(full)
         if (i /= 0) then
            s = to_lower( trim(adjustl(self%store_strs(i)%val)) )
            if (s=='true' .or. s=='false') then
               tf = .true.
            else if (rel) then
               if (s=='yes' .or. s=='no' .or. s=='on' .or. s=='off' .or. s=='1' .or. s=='0') tf = .true.
            end if
         end if
      end block
   end function cfg_is_boolean


   logical function cfg_is_string(self, key, include_empty_scalar) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      logical, optional,   intent(in) :: include_empty_scalar  ! default .false. here
      character(len=PARAMLEN) :: full
      integer :: i
      logical :: incl

      tf = .false.
      incl = .false.; if (present(include_empty_scalar)) incl = include_empty_scalar

      full = self%resolve_key(key)
      if (len_trim(full) == 0) return
      if (self%is_null(key, include_empty_scalar=incl)) return

      i = self%find_str_key(full)
      if (i /= 0) then
         ! It's a stored string (including empty "" if incl=.true.)
         tf = .true.; return
      end if
      ! Anything else (real/int/bool/map/seq) → not string
   end function cfg_is_string

   logical function cfg_is_disabled(self, key) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),    intent(in) :: key
      character(len=PARAMLEN) :: full
      integer :: i
      character(len=STRLEN) :: s

      ! Null / empty / "null"/"none"/"~" → disabled
      tf = self%is_null(key, include_empty_scalar=.true., coerce_none=.true.)
      if (tf) return

      full = self%resolve_key(key)
      if (len_trim(full) == 0) then
         tf = .true.; return  ! unresolved → treat as disabled
      end if

      ! Logical false → disabled
      i = self%find_log_key(full)
      if (i /= 0) then
         tf = .not. self%store_logs(i)%val
         return
      end if

      ! String false-tokens → disabled (case-insensitive)
      i = self%find_str_key(full)
      if (i /= 0) then
         s = to_lower(trim(adjustl(self%store_strs(i)%val)))
         if (s == 'off' .or. s == 'no' .or. s == 'false') then
            tf = .true.; return
         end if
      end if

      ! Numeric nodes (including 0) → NOT disabled
      tf = .false.
   end function cfg_is_disabled




   !==================== Resolution & lookups ====================

   function resolve_key(self, query) result(fullkey)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: query
      character(len=PARAMLEN) :: fullkey
      integer :: i, hits
      character(len=PARAMLEN) :: q

      q = trim(query)
      fullkey = ''

      if (self%exists_fullkey(q)) then
         fullkey = q; return
      end if

      hits = 0
      if (allocated(self%store_reals)) then
         do i=1,size(self%store_reals)
            if (basename(self%store_reals(i)%key) == q) then
               hits = hits + 1; fullkey = self%store_reals(i)%key
            end if
         end do
      end if
      if (allocated(self%store_ints)) then
         do i=1,size(self%store_ints)
            if (basename(self%store_ints(i)%key) == q) then
               hits = hits + 1; fullkey = self%store_ints(i)%key
            end if
         end do
      end if
      if (allocated(self%store_logs)) then
         do i=1,size(self%store_logs)
            if (basename(self%store_logs(i)%key) == q) then
               hits = hits + 1; fullkey = self%store_logs(i)%key
            end if
         end do
      end if
      if (allocated(self%store_strs)) then
         do i=1,size(self%store_strs)
            if (basename(self%store_strs(i)%key) == q) then
               hits = hits + 1; fullkey = self%store_strs(i)%key
            end if
         end do
      end if

      if (allocated(self%store_nulls)) then
         do i=1,size(self%store_nulls)
            if (basename(self%store_nulls(i)%key) == q) then
               hits = hits + 1; fullkey = self%store_nulls(i)%key
            end if
         end do
      end if


      if (hits > 1) call stop_missingpar('resolve_key','Parameter name "'//q//'" is ambiguous; use full dotted key.')
   end function resolve_key

   function exists_fullkey(self, k) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      logical :: tf
      tf = (self%find_real_key(k) /= 0) .or. &
           (self%find_int_key(k)  /= 0) .or. &
           (self%find_log_key(k)  /= 0) .or. &
           (self%find_str_key(k)  /= 0) .or. &
           (self%find_null_key(k) /= 0 )
   end function exists_fullkey

   function lookup_string(self, k) result(val)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      character(len=STRLEN) :: val
      integer :: i
      val = ''
      i = self%find_str_key(k)
      if (i /= 0) val = self%store_strs(i)%val
   end function lookup_string

   !==================== Append & find helpers ====================

   subroutine append_real(self, k, v)
      class(ConfigParams), intent(inout) :: self
      character(len=*),   intent(in)    :: k
      real(rk),           intent(in)    :: v
      type(cfg_real_t), allocatable :: tmp(:)
      integer :: n
      n = 0; if (allocated(self%store_reals)) n = size(self%store_reals)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = self%store_reals
      tmp(n+1)%key = left_trim_to(k, PARAMLEN)
      tmp(n+1)%val = v
      call move_alloc(tmp, self%store_reals)
   end subroutine append_real

   subroutine append_int(self, k, v)
      class(ConfigParams), intent(inout) :: self
      character(len=*),   intent(in)    :: k
      integer,            intent(in)    :: v
      type(cfg_int_t), allocatable :: tmp(:)
      integer :: n
      n = 0; if (allocated(self%store_ints)) n = size(self%store_ints)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = self%store_ints
      tmp(n+1)%key = left_trim_to(k, PARAMLEN)
      tmp(n+1)%val = v
      call move_alloc(tmp, self%store_ints)
   end subroutine append_int

   subroutine append_log(self, k, v)
      class(ConfigParams), intent(inout) :: self
      character(len=*),   intent(in)    :: k
      logical,            intent(in)    :: v
      type(cfg_logical_t), allocatable :: tmp(:)
      integer :: n
      n = 0; if (allocated(self%store_logs)) n = size(self%store_logs)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = self%store_logs
      tmp(n+1)%key = left_trim_to(k, PARAMLEN)
      tmp(n+1)%val = v
      call move_alloc(tmp, self%store_logs)
   end subroutine append_log

   subroutine append_str(self, k, v)
      class(ConfigParams), intent(inout) :: self
      character(len=*),   intent(in)    :: k
      character(len=*),   intent(in)    :: v
      type(cfg_string_t), allocatable :: tmp(:)
      integer :: n
      n = 0; if (allocated(self%store_strs)) n = size(self%store_strs)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = self%store_strs
      tmp(n+1)%key = left_trim_to(k, PARAMLEN)
      tmp(n+1)%val = left_trim_to(v, STRLEN)
      call move_alloc(tmp, self%store_strs)
   end subroutine append_str

   function find_real_key(self, k) result(idx)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      integer :: idx, i
      idx = 0
      if (.not.allocated(self%store_reals)) return
      do i=1, size(self%store_reals)
         if (trim(self%store_reals(i)%key) == trim(k)) then
            idx = i; return
         end if
      end do
   end function find_real_key

   function find_int_key(self, k) result(idx)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      integer :: idx, i
      idx = 0
      if (.not.allocated(self%store_ints)) return
      do i=1, size(self%store_ints)
         if (trim(self%store_ints(i)%key) == trim(k)) then
            idx = i; return
         end if
      end do
   end function find_int_key

   function find_log_key(self, k) result(idx)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      integer :: idx, i
      idx = 0
      if (.not.allocated(self%store_logs)) return
      do i=1, size(self%store_logs)
         if (trim(self%store_logs(i)%key) == trim(k)) then
            idx = i; return
         end if
      end do
   end function find_log_key

   function find_str_key(self, k) result(idx)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      integer :: idx, i
      idx = 0
      if (.not.allocated(self%store_strs)) return
      do i=1, size(self%store_strs)
         if (trim(self%store_strs(i)%key) == trim(k)) then
            idx = i; return
         end if
      end do
   end function find_str_key

   subroutine append_null(self, k)
      class(ConfigParams), intent(inout) :: self
      character(len=*),    intent(in)    :: k
      type(cfg_null_t), allocatable :: tmp(:)
      integer :: n
      n = 0; if (allocated(self%store_nulls)) n = size(self%store_nulls)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = self%store_nulls
      tmp(n+1)%key = left_trim_to(k, PARAMLEN)
      call move_alloc(tmp, self%store_nulls)
   end subroutine append_null

   integer function find_null_key(self, k) result(idx)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      integer :: i
      idx = 0
      if (.not. allocated(self%store_nulls)) return
      do i=1,size(self%store_nulls)
         if (self%store_nulls(i)%key == k) then; idx = i; return; end if
      end do
   end function find_null_key

   !==================== Utils ====================

   pure logical function startswith(s, prefix)
      character(len=*), intent(in) :: s, prefix
      integer :: lp
      lp = len_trim(prefix)
      if (lp == 0) then
         startswith = .true.
      else if (len_trim(s) < lp) then
         startswith = .false.
      else
         startswith = (s(1:lp) == prefix(1:lp))
      end if
   end function startswith

   pure function basename(k) result(b)
      character(len=*), intent(in) :: k
      character(len=PARAMLEN) :: b
      integer :: i, lt
      b = ''; lt = len_trim(k); i = lt
      do while (i >= 1)
         if (k(i:i) == '.') exit
         i = i - 1
      end do
      if (i < 1) then
         b(1:lt) = k(1:lt)
      else
         if (i < lt) b(1:lt-i) = k(i+1:lt)
      end if
   end function basename

   pure logical function is_true_string(s)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      t = to_lower(trim(s))
      is_true_string = (t=='true' .or. t=='t' .or. t=='yes' .or. t=='y' .or. t=='on' .or. t=='1')
   end function is_true_string

   pure logical function is_false_string(s)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      t = to_lower(trim(s))
      is_false_string = (t=='false' .or. t=='f' .or. t=='no' .or. t=='n' .or. t=='off' .or. t=='0')
   end function is_false_string

   pure function left_trim_to(s, maxlen) result(out)
      character(len=*), intent(in) :: s
      integer, intent(in) :: maxlen
      character(len=maxlen) :: out
      integer :: lt, n
      out = ''
      lt = len_trim(s)
      n = min(maxlen, lt)
      if (n > 0) out(1:n) = s(1:n)
   end function left_trim_to

   subroutine stop_missingpar(where, msg)
      character(len=*), intent(in) :: where, msg
      write(*,'(A,1X,A)') '[FATAL:'//trim(where)//']', trim(msg)
      stop 1
   end subroutine stop_missingpar

   pure function limit_msg_num(key, op, limit, val) result(msg)
      character(*), intent(in) :: key, op
      real(rk),    intent(in) :: limit, val
      character(:), allocatable :: msg
      character(48) :: slimit, sval
      write(slimit,'(G0)') limit
      write(sval,  '(G0)') val
      msg = '"'//trim(key)//'" ('//trim(sval)//') must be '//trim(op)//' '//trim(slimit)//'.'
   end function limit_msg_num



   pure function join_key(prefix, key) result(out)
      character(len=*), intent(in) :: prefix, key
      character(len=PARAMLEN) :: out
      out = ''
      if (len_trim(prefix) == 0) then
         out = left_trim_to(key, PARAMLEN)
      else
         out = left_trim_to(trim(prefix)//'.'//trim(key), PARAMLEN)
      end if
   end function join_key


end module read_config_yaml
