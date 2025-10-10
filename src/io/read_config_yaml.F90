module read_config_yaml
   !! Reads and loads parameters from a YAML config file using FABM's yaml parsing libraries.
   !! - Supports full dotted parameters ("physics.time.dt") and unique name-only parameters ("dt")
   !! - Numeric: returns real(rk); parses ints/reals/boolean-like strings
   !! - Logical: true/false/yes/no/on/off/0/1 (case-insensitive)
   !!
   !! Usage:
   !!   use read_config_yaml
   !!   type(ConfigParams) :: params
   !!   call params%init()
   !!   call params%load_yaml_content('model.yaml',    'model')
   !!   call params%load_yaml_content('physics.yaml',  'physics')
   !!   call params%load_yaml_content('transport.yaml','transport')
   !!   params%get_param_num('physics.level1.param',default_value)
   !!
   !!   real(rk) :: dt
   !!   dt = params%get_param_num('physics.time.dt')
   !!   print *, params%get_param_str('output.file', 'out.nc')   ! name-only if unique
   !!
   use precision_types, only: rk
   use yaml_types
   use yaml, only: yaml_parse => parse, yaml_error_length => error_length

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

   type :: ConfigParams
      private
      type(cfg_real_t),    allocatable :: store_reals(:)
      type(cfg_int_t),     allocatable :: store_ints(:)
      type(cfg_logical_t), allocatable :: store_logs(:)
      type(cfg_string_t),  allocatable :: store_strs(:)
   contains
      procedure :: init                => cfg_init
      procedure :: clear               => cfg_clear
      procedure :: load_yaml_content   => cfg_load_yaml_content
      procedure :: has_param           => cfg_has_param
      procedure :: list_params         => cfg_list_params
      procedure :: get_param_str       => cfg_get_param_str
      procedure :: get_param_num       => cfg_get_param_num
      procedure :: get_param_logical   => cfg_get_param_logical
      ! internals
      procedure, private :: walk_dictionary
      procedure, private :: store_scalar
      procedure, private :: resolve_key
      procedure, private :: exists_exact
      procedure, private :: lookup_string
      procedure, private :: append_real
      procedure, private :: append_int
      procedure, private :: append_log
      procedure, private :: append_str
      procedure, private :: find_real_key
      procedure, private :: find_int_key
      procedure, private :: find_log_key
      procedure, private :: find_str_key
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
            call stop_missingpar('walk_dictionary', trim(node%path)//' must not be null.')
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

      ! 1) Try real (covers ints in many builds too)
      rval = s%to_real(default=0.0_rk, success=ok)
      if (ok) then
         call self%append_real(flat_key, rval); return
      end if

      ! 2) Try parse integer from string
      str = adjustl(trim(s%string))
      read(str, *, iostat=ios) ival
      if (ios == 0) then
         call self%append_int(flat_key, ival); return
      end if

      ! 3) Try boolean-like strings
      low = to_lower(str)
      if (is_true_string(low))  then; call self%append_log(flat_key, .true. ); return; end if
      if (is_false_string(low)) then; call self%append_log(flat_key, .false.); return; end if

      ! 4) Fallback: store as string
      call self%append_str(flat_key, str)
   end subroutine store_scalar

   !==================== Getters ====================

   function cfg_get_param_str(self, key, default) result(val)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      character(len=*),   intent(in), optional :: default
      character(len=STRLEN) :: val
      character(len=PARAMLEN) :: full

      full = self%resolve_key(key)
      if (len_trim(full) /= 0) then
         val = self%lookup_string(full)
      else
         if (present(default)) then
            val = default
         else
            call stop_missingpar('get_param_str', 'Missing string parameter: '//trim(key))
         end if
      end if
   end function cfg_get_param_str

   function cfg_get_param_logical(self, key, default) result(out)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      logical,            intent(in), optional :: default
      logical :: out
      character(len=PARAMLEN) :: full
      character(len=STRLEN) :: s
      integer :: i

      full = self%resolve_key(key)
      if (len_trim(full) == 0) then
         if (present(default)) then
            out = default
         else
            call stop_missingpar('get_param_logical', 'Missing logical parameter: '//trim(key))
         end if
         return
      end if

      i = self%find_log_key(full); if (i /= 0) then; out = self%store_logs(i)%val; return; end if
      i = self%find_int_key(full); if (i /= 0) then; out = (self%store_ints(i)%val /= 0); return; end if
      i = self%find_real_key(full);if (i /= 0) then; out = (abs(self%store_reals(i)%val) > 0.5_rk); return; end if

      s = self%lookup_string(full)
      s = to_lower(trim(s))
      if (is_true_string(s)) then
         out = .true.
      else if (is_false_string(s)) then
         out = .false.
      else
         out = (len_trim(s) > 0)   ! last-resort heuristic
      end if
   end function cfg_get_param_logical

   function cfg_get_param_num(self, key, default) result(val)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      real(rk),           intent(in), optional :: default
      real(rk) :: val
      character(len=PARAMLEN) :: full
      character(len=STRLEN) :: s
      integer :: i, ios, itmp
      real(rk) :: rtmp

      full = self%resolve_key(key)
      if (len_trim(full) == 0) then
         if (present(default)) then
            val = default
         else
            call stop_missingpar('get_param_num', 'Missing numeric parameter: '//trim(key))
         end if
         return
      end if

      i = self%find_real_key(full); if (i /= 0) then; val = self%store_reals(i)%val; return; end if
      i = self%find_int_key(full);  if (i /= 0) then; val = real(self%store_ints(i)%val, rk); return; end if
      i = self%find_log_key(full);  if (i /= 0) then; val = merge(1.0_rk, 0.0_rk, self%store_logs(i)%val); return; end if

      s = adjustl(trim(self%lookup_string(full)))
      read(s, *, iostat=ios) itmp
      if (ios == 0) then; val = real(itmp, rk); return; end if
      read(s, *, iostat=ios) rtmp
      if (ios == 0) then; val = rtmp; return; end if

      call stop_missingpar('get_param_num', 'Parameter "'//trim(key)//'" is not numeric.')
   end function cfg_get_param_num

   !==================== Queries ====================

   function cfg_has_param(self, key) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: key
      logical :: tf
      character(len=PARAMLEN) :: full
      full = self%resolve_key(key)
      tf = (len_trim(full) /= 0)
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

      allocate(list(max(1,n)))
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
   end function cfg_list_params

   !==================== Resolution & lookups ====================

   function resolve_key(self, query) result(fullkey)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: query
      character(len=PARAMLEN) :: fullkey
      integer :: i, hits
      character(len=PARAMLEN) :: q

      q = trim(query)
      fullkey = ''

      if (self%exists_exact(q)) then
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

      if (hits > 1) call stop_missingpar('resolve_key','Parameter name "'//q//'" is ambiguous; use full dotted key.')
   end function resolve_key

   function exists_exact(self, k) result(tf)
      class(ConfigParams), intent(in) :: self
      character(len=*),   intent(in) :: k
      logical :: tf
      tf = (self%find_real_key(k) /= 0) .or. (self%find_int_key(k) /= 0) .or. &
           (self%find_log_key(k)  /= 0) .or. (self%find_str_key(k)  /= 0)
   end function exists_exact

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

   pure function to_lower(s) result(o)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: o
      integer :: i, ia
      do i=1,len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('A') .and. ia <= iachar('Z')) then
            o(i:i) = achar(ia + 32)
         else
            o(i:i) = s(i:i)
         end if
      end do
   end function to_lower

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
