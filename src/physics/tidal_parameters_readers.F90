module tidal_parameters_readers
  use iso_fortran_env,  only: real64
  use read_config_yaml, only: ConfigParams
  implicit none
  private
  public :: read_tidal_parameters

  integer, parameter :: rk = selected_real_kind(12, 200)

  !-------------DATA STRUCTURES------------------------------------------
  !Data structure for each Tidal constituent
  type, public :: Constituent
    character(len=:), allocatable :: name
    real(real64) :: sema    = 0.0_real64
    real(real64) :: semi    = 0.0_real64
    real(real64) :: inc_deg = 0.0_real64
    real(real64) :: pha_deg = 0.0_real64
  end type
  ! Data structure for all Tidal Parameters
  type, public :: TidalParams
    type(Constituent), allocatable :: c(:)  ! all constituents
  contains
    procedure :: get => get_constituent
  end type TidalParams
  

contains

  !--------------------- Main subroutine------------------------------------
  subroutine read_tidal_parameters(cfg, tide_const)
    !! Reads tidal ellipses either from file (preferred) or from YAML block.
    type(ConfigParams), intent(in)  :: cfg
    type(TidalParams),  intent(out) :: tide_const

    logical :: from_file
    from_file = cfg%get_param_logical('physics.tides.read_from_file', .true.)

    if (from_file) then
      call read_from_file(cfg, tide_const)
    else
      call read_from_yaml(cfg, tide_const)
    end if
  end subroutine read_tidal_parameters

  !------------------ Read tidal parameters from dat file ----------------------
  subroutine read_from_file(cfg, tide_const)
    type(ConfigParams), intent(in)  :: cfg
    type(TidalParams),  intent(out) :: tide_const

    character(len=:), allocatable :: fname
    real(real64) :: lon, lat, tol
    integer :: idx, n, k

    ! temps to receive data from the file reader
    character(len=:), allocatable :: names(:)
    real(real64),     allocatable :: sema(:), semi(:), inc(:), pha(:)

    fname = cfg%get_param_str('physics.tides.filename')
    tol   = cfg%get_param_num('physics.tides.tol_deg', 0.1_real64)
    lat   = cfg%get_param_num('main.location.latitude')
    lon   = cfg%get_param_num('main.location.longitude')

    call read_tides_at_point(fname, lon, lat, names, sema, semi, inc, pha, idx, tol_deg=tol)

    ! Fill each Constituent container
    n = size(names)
    if (n <= 0) stop 'read_from_file: no constituents returned.'

    allocate(tide_const%c(n))
    do k = 1, n
      tide_const%c(k)%name    = trim(names(k))
      tide_const%c(k)%sema    = sema(k)
      tide_const%c(k)%semi    = semi(k)
      tide_const%c(k)%inc_deg = inc(k)
      tide_const%c(k)%pha_deg = pha(k)
    end do
  end subroutine read_from_file

  !---------------- Read tidal parameters from the yaml file ----------------------
  subroutine read_from_yaml(cfg, tide_const)
    implicit none
    type(ConfigParams), intent(in)  :: cfg
    type(TidalParams),  intent(out) :: tide_const

    ! Supported list (extend as needed)
    character(len=*), parameter :: constituents(5) = ['M2','S2','K1','O1','N2']

    integer :: i, nall
    character(len=:), allocatable :: base
    real(real64), parameter :: MISS = -huge(1.0_real64)
    real(real64) :: a, b, ai, ap
    logical :: found

    nall = size(constituents)
    allocate(tide_const%c(nall))

    do i = 1, nall
      base = 'physics.tides.ellipse_params.' // trim(constituents(i))

      ! Retrieve with sentinel defaults
      a  = cfg%get_param_num(base // '.semi_major', MISS)
      b  = cfg%get_param_num(base // '.semi_minor', MISS)
      ai = cfg%get_param_num(base // '.inclination', 0.0_real64)
      ap = cfg%get_param_num(base // '.phase_ang',  0.0_real64)

      found = (a /= MISS .and. b /= MISS)

      if (found) then
        tide_const%c(i)%name    = lower(trim(constituents(i)))
        tide_const%c(i)%sema    = a
        tide_const%c(i)%semi    = b
        tide_const%c(i)%inc_deg = ai
        tide_const%c(i)%pha_deg = ap
      else
        ! Constituent missing -> fill zeros, warn
        tide_const%c(i)%name    = lower(trim(constituents(i)))
        tide_const%c(i)%sema    = 0.0_real64
        tide_const%c(i)%semi    = 0.0_real64
        tide_const%c(i)%inc_deg = 0.0_real64
        tide_const%c(i)%pha_deg = 0.0_real64
        write(*,'(A)') 'Warning [read_from_yaml]: Constituent ' // trim(constituents(i)) // &
                      ' not found in YAML; set to zero.'
      end if
    end do
  end subroutine read_from_yaml 
  

  !================HELPERS====================================================

  pure function lower(s) result(t)
    ! Returns the lowercase version of the string s
    ! Converts character by character
    character(len=*), intent(in) :: s
    character(len=len(s))        :: t
    integer :: i
    t = s
    do i=1,len(s)
      select case(s(i:i))
      case('A':'Z'); t(i:i) = achar(iachar(s(i:i)) + 32)
      case default
      end select
    end do
  end function lower

  subroutine get_column_names(line, col_names)
    !! Split a header line on whitespace (compressing multiple spaces)
    !! Return an array with column names
    character(len=*), intent(in) :: line
    character(len=:), allocatable, intent(out) :: col_names(:)

    character(len=:), allocatable :: s
    integer :: L, i, istart, n, maxlen

    s = adjustl(line(1:len_trim(line)))
    L = len_trim(s)
    istart = 0
    n = 0
    maxlen = 0

    ! Count found names and find max name length
    do i = 1, L+1
      if (i==L+1 .or. s(i:i)==' ') then
        if (istart>0) then
          n = n + 1
          maxlen = max(maxlen, i - istart)
          istart = 0
        end if
      else
        if (istart==0) istart = i
      end if
    end do

    allocate(character(len=maxlen) :: col_names(n))

    ! Copy column names
    istart = 0
    n = 0
    do i = 1, L+1
      if (i==L+1 .or. s(i:i)==' ') then
        if (istart>0) then
          n = n + 1
          col_names(n) = ' '
          col_names(n)(1:i-istart) = s(istart:i-1)
          istart = 0
        end if
      else
        if (istart==0) istart = i
      end if
    end do
  end subroutine get_column_names

  subroutine parse_header(unit, nconst, names, idx0)
    !! Read a header line starting with '#'.
    !! Outputs:
    !!  nconst          : number of constituents (groups of four after lon/lat)
    !!  names(1:nconst) : constituent names (lowercase, e.g. 'm2','k1')
    !!  idx0(1:nconst)  : 1-based starting column index for SEMA in numeric rows
    integer, intent(in) :: unit
    integer, intent(out) :: nconst
    character(len=:), allocatable, intent(out) :: names(:)
    integer, allocatable, intent(out) :: idx0(:)

    character(len=2048) :: line
    character(len=:), allocatable :: cols(:)
    integer :: ncol, k, ius, max_name_len

    rewind(unit)
    do
      read(unit,'(A)',end=900) line
      if (len_trim(line)==0) cycle
      if (line(1:1) == '#') exit
    end do

    call get_column_names(trim(line(2:)), cols)  ! drop the leading '#'
    ncol = size(cols)

    if (ncol < 2) stop 'parse_header: header must include "lon lat".'
    if (lower(cols(1)) /= 'lon' .or. lower(cols(2)) /= 'lat') then
      stop 'parse_header: first tokens must be "lon lat".'
    end if
    if (mod(ncol-2, 4) /= 0) then
      stop 'parse_header: columns after lon/lat must be in groups of 4.'
    end if

    nconst = (ncol - 2) / 4
    allocate(idx0(nconst))

    ! Determine max constituent-name length to allocate once
    max_name_len = 0
    do k = 1, nconst
      ius = index(cols(2 + 4*(k-1) + 1), '_')       ! "<name>_SEMA"
      if (ius <= 1) stop 'parse_header: malformed token; expected "<name>_SEMA".'
      max_name_len = max(max_name_len, ius - 1)
    end do
    allocate(character(len=max_name_len) :: names(nconst))

    ! Fill names and starting indices (SEMA column)
    do k = 1, nconst
      ius = index(cols(2 + 4*(k-1) + 1), '_')
      names(k) = ' '
      names(k)(1:ius-1) = lower(adjustl(cols(2 + 4*(k-1) + 1)(1:ius-1)))
      idx0(k)  = 2 + 1 + 4*(k-1)    ! 1-based numeric column index for SEMA
    end do
    return
  900 continue
    stop 'parse_header: header line starting with "#" not found.'
  end subroutine parse_header

  subroutine read_tides_at_point(filename, lon_t, lat_t, &
                                   names, sema, semi, inc, pha, found_idx, tol_deg)
    !! Extract the tidal parameters corresponding to the closest gridpoint                              
    !! Read file, parse header, select nearest or exact lon/lat row,
    !! return arrays in header order.
    character(len=*), intent(in) :: filename
    real(rk), intent(in) :: lon_t, lat_t
    character(len=:), allocatable, intent(out) :: names(:)
    real(rk), allocatable, intent(out) :: sema(:), semi(:), inc(:), pha(:)
    integer, intent(out) :: found_idx        ! 1-based data row index
    real(rk), intent(in), optional :: tol_deg
    real(rk), parameter :: eps_exact = 5.0e-4_rk

    integer :: iu, io, nconst, ncols, i, k
    character(len=4096) :: line
    integer, allocatable :: idx0(:)
    real(rk), allocatable :: vals(:), best_vals(:)
    real(rk) :: lon, lat, d2, best_d2, tol
    integer :: data_row_count

    open(newunit=iu, file=filename, status='old', action='read', iostat=io)
    if (io /= 0) then
        write(*,'(A,1X,A,1X,A)') 'File: ', trim(filename), 'with tidal parameters not found or not readable.'
        stop
    end if

    call parse_header(iu, nconst, names, idx0)
    ncols = 2 + 4*nconst
    allocate(vals(ncols))

    tol = 5.0e-4_rk       ; if (present(tol_deg))    tol = tol_deg
    

    best_d2 = huge(1.0_rk)
    found_idx = 0
    data_row_count = 0
    if (allocated(best_vals)) deallocate(best_vals)

    do
      read(iu,'(A)', end=100) line
      if (len_trim(line)==0) cycle
      if (line(1:1) == '#')  cycle

      vals = 0.0_rk
      read(line,*,iostat=io) (vals(i), i=1, ncols)
      if (io /= 0) cycle

      lon = vals(1); lat = vals(2)
      data_row_count = data_row_count + 1

      if (tol == 0.0_rk) then
        ! "exact" match to printed precision
        if (abs(lon - lon_t) <= eps_exact .and. abs(lat - lat_t) <= eps_exact) then
          found_idx = data_row_count
          allocate(best_vals(ncols))
          best_vals = vals
          exit
        end if
      else
        ! nearest neighbor (Euclidean in degrees)
        d2 = (lon - lon_t)**2 + (lat - lat_t)**2
        if (d2 < best_d2) then
          best_d2 = d2
          if (.not. allocated(best_vals)) allocate(best_vals(ncols))
          best_vals = vals
          found_idx = data_row_count
        end if
      end if
    end do
100 continue
    close(iu)

    if (tol == 0.0_rk) then
      if (found_idx == 0) stop 'No exact tidal-parameters match at this lon/lat (to 3 decimal points).'
    else
      if (found_idx == 0) stop 'No data rows found.'
      if (sqrt(best_d2) > tol) stop 'Nearest point for tidal parameters is outside tol_deg.'
    end if

    allocate(sema(nconst), semi(nconst), inc(nconst), pha(nconst))
    do k = 1, nconst
      sema(k) = best_vals(idx0(k)    )
      semi(k) = best_vals(idx0(k) + 1)
      inc(k)  = best_vals(idx0(k) + 2)
      pha(k)  = best_vals(idx0(k) + 3)
    end do
  end subroutine read_tides_at_point  

  !----------- TidalParams Method to access each constituent values
  pure subroutine get_constituent(self, name, c)
    class(TidalParams), intent(in)  :: self
    character(len=*),    intent(in) :: name
    type(Constituent),   intent(out):: c
    integer :: i
    character(len=:), allocatable :: key

    key = lower(trim(name))
    c = Constituent()  

    do i = 1, size(self%c)
      if (trim(self%c(i)%name) == key) then
        c = self%c(i)
        return
      end if
    end do
  end subroutine get_constituent

end module tidal_parameters_readers
