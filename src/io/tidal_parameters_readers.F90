module tidal_readers
  use precision_types,  only: rk
  use read_config_yaml, only: ConfigParams
  use precision_utils,  only: is_equal, is_unequal, is_zero
  use tidal,            only: TidalConstituent, TidalSet
  use str_utils,        only: to_lower, to_upper
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: read_tidal_parameters  

contains

  !--------------------- Main subroutine------------------------------------
  subroutine read_tidal_parameters(cfg,lat,lon,tide_set)
    !! Reads tidal ellipses either from file (preferred) or from YAML block.
    type(ConfigParams), intent(in)     :: cfg
    real(rk),           intent(in)     :: lat,lon
    type(TidalSet),     intent(out)    :: tide_set
    


    logical :: from_file
    from_file = cfg%get_param_logical('tides.read_from_file', default=.true., strict=.false.)

    if (from_file) then
      call read_from_file(cfg,lat,lon, tide_set)
    else
      call read_from_yaml(cfg, tide_set)
    end if
  end subroutine read_tidal_parameters

  !------------------ Read tidal parameters from dat file ----------------------
  subroutine read_from_file(cfg,lat,lon, tide_set)
    type(ConfigParams), intent(in)  :: cfg
    real(rk),           intent(in)  :: lat,lon
    type(TidalSet),  intent(out) :: tide_set    

    character(len=:), allocatable :: fname
    real(rk) :: tol
    integer  :: idx, n, k

    ! temps to receive data from the file reader
    character(len=:), allocatable :: names(:)
    real(rk),     allocatable :: smaj(:), smin(:), inc(:), pha(:)

    fname = cfg%get_param_str('tides.filename', required=.true., trim_value=.true., empty_ok=.false.)
    tol   = cfg%get_param_num('tides.tol_deg', 0.1_rk, nonnegative=.true., max=10.0_rk, finite=.true.)

    call read_tides_at_point(fname, lon, lat, names, smaj, smin, inc, pha, idx, tol_deg=tol)

    ! Fill each Constituent container
    n = size(names)
    if (n <= 0) stop 'read_from_file: no constituents returned.'

    allocate(tide_set%C(n))
    do k = 1, n
      tide_set%C(k)%name    = to_upper(trim(names(k)))
      tide_set%C(k)%smaj    = smaj(k)
      tide_set%C(k)%smin    = smin(k)
      tide_set%C(k)%theta   = inc(k)
      tide_set%C(k)%phase   = pha(k)
    end do
  end subroutine read_from_file

  !---------------- Read tidal parameters from the yaml file ----------------------
  subroutine read_from_yaml(cfg, tide_set)
    implicit none
    type(ConfigParams), intent(in)  :: cfg
    type(TidalSet),  intent(out) :: tide_set

    ! Supported list (extend as needed)
    character(len=*), parameter :: constituents(5) = ['M2','S2','K1','O1','N2']

    integer :: i, nall
    character(len=:), allocatable :: base
    real(rk), parameter :: MISS = -huge(1.0_rk)
    real(rk) :: a, b, ai, ap
    logical  :: fa, fb, fai, fap, found   ! found flags

    nall = size(constituents)
    allocate(tide_set%c(nall))

    do i = 1, nall
      base = 'tides.ellipse_params.' // trim(constituents(i))

      ! Retrieve cnstituents
      a  = cfg%get_param_num(base // '.semi_major', found=fa, nonnegative=.true., finite=.true.)
      b  = cfg%get_param_num(base // '.semi_minor', found=fb, finite=.true.)
      ai = cfg%get_param_num(base // '.inclination', default=0.0_rk, found=fai, finite=.true.)
      ap = cfg%get_param_num(base // '.phase_ang', default=0.0_rk, found=fap, finite=.true.)

      found = fa .and. fb

      if (found) then
        if (abs(b) > a + 1.0e-12_rk) then
          error stop 'Tidal ellipse invalid at '//trim(base)//': |smin_minor| > smin_major.'
        end if

        tide_set%c(i)%name    = to_upper(trim(constituents(i)))
        tide_set%c(i)%smaj    = a
        tide_set%c(i)%smin    = b
        tide_set%c(i)%theta = ai
        tide_set%c(i)%phase = ap
      else
        ! Constituent missing -> fill zeros, warn
        tide_set%c(i)%name    = to_upper(trim(constituents(i)))
        tide_set%c(i)%smaj    = 0.0_rk
        tide_set%c(i)%smin    = 0.0_rk
        tide_set%c(i)%theta = 0.0_rk
        tide_set%c(i)%phase = 0.0_rk
        write(*,'(A)') 'Warning [read_from_yaml]: Constituent ' // trim(constituents(i)) // &
                      ' not found in YAML; set to zero.'
      end if
    end do
  end subroutine read_from_yaml 
  

  !================HELPERS====================================================

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
    !!  names(1:nconst) : constituent names 
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
    if (to_lower(cols(1)) /= 'lon' .or. to_lower(cols(2)) /= 'lat') then
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
      names(k)(1:ius-1) = to_upper(adjustl(cols(2 + 4*(k-1) + 1)(1:ius-1)))
      idx0(k)  = 2 + 1 + 4*(k-1)    ! 1-based numeric column index for SEMA
    end do
    return
  900 continue
    stop 'parse_header: header line starting with "#" not found.'
  end subroutine parse_header

  subroutine read_tides_at_point(filename, lon_t, lat_t, &
                                   names, smaj, smin, inc, pha, found_idx, tol_deg)
    !! Extract the tidal parameters corresponding to the closest gridpoint                              
    !! Read file, parse header, select nearest or exact lon/lat row,
    !! return arrays in header order.
    character(len=*), intent(in) :: filename
    real(rk), intent(in) :: lon_t, lat_t
    character(len=:), allocatable, intent(out) :: names(:)
    real(rk), allocatable, intent(out) :: smaj(:), smin(:), inc(:), pha(:)
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

      if (is_zero(tol, abs_=5.0e-4_rk)) then
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

    if (is_zero(tol, abs_=5.0e-4_rk)) then
      if (found_idx == 0) stop 'No exact tidal-parameters match at this lon/lat (to 3 decimal points).'
    else
      if (found_idx == 0) stop 'No data rows found.'
      if (sqrt(best_d2) > tol) stop 'Nearest point for tidal parameters is outside tol_deg.'
    end if

    allocate(smaj(nconst), smin(nconst), inc(nconst), pha(nconst))
    do k = 1, nconst
      smaj(k) = best_vals(idx0(k)    )
      smin(k) = best_vals(idx0(k) + 1)
      inc(k)  = best_vals(idx0(k) + 2)
      pha(k)  = best_vals(idx0(k) + 3)
    end do
  end subroutine read_tides_at_point   

end module tidal_readers
