!=======================================================================================
!   Shelf Sea Physics:
!  1-D MODEL OF THE EQUATION OF MOTION USING THE Canuto k-e TURBULENCE CLOSURE SCHEME
!
!   This code is fully based on the S2P3 model lineage:
!
!   • Original S2P3 (v7.0):
!       Jonathan Sharples – Univ. of Liverpool & NOC
!       Documented in: Simpson & Sharples (2012), CUP.
!
!   • Regional framework S2P3-R (v1.0):
!       Marsh, Hickman & Sharples (2015), GMD 8, 3163–3178.
!
!   • Large-scale/efficient S2P3-R v2.0:
!       Halloran et al. (2021), GMD 14, 6177–6195.
!
!=======================================================================================
module tidal
  use array_utils,      only: str_to_real_vec
  use find_utils,       only: argmin_abs_vec
  use geo_utils,        only: simple_distance_deg
  use precision_types,  only: rk
  use physics_params,   only: Omega  
  use read_config_yaml, only: ConfigParams
  use str_utils,        only: to_upper
  use text_parser,      only: TextTable, parse_text_table
  use trigonometrics,   only: deg2rad
  
  implicit none
  private

  public :: TidalConstituent, TidalSet 
  
  type :: TidalConstituent
     character(len=8) :: name = ''     
     real(rk) :: smaj  = 0.0_rk      ! m s^-1
     real(rk) :: smin  = 0.0_rk      ! m s^-1 
     real(rk) :: theta = 0.0_rk      ! rad (inclination of semi-major, CCW from East) -> Orient in S2P3
     real(rk) :: phase = 0.0_rk      ! rad   -> major_phase in S2P3
     real(rk) :: omega = 0.0_rk      ! rad s^-1
     real(rk) :: polarisn = 0.0_rk   ! polarization (ellipse axis ratio)
     real(rk) :: slopex = 0.0_rk     ! constituent acceleration amplitudes along the ellipse major axis (m s-2)
     real(rk) :: slopey = 0.0_rk     ! constituent acceleration amplitudes along the ellipse minor axis (m s-2)
     real(rk) :: ctheta = 1.0_rk     ! Cos of theta for the rotation matrix (calculate once and avoid calculating them each step later)
     real(rk) :: stheta = 0.0_rk     ! Sin of theta
  end type

  type :: TidalSet
     type(TidalConstituent), allocatable :: C(:)
     real(rk) :: f0      = 0.0_rk     ! 2*Omega*sin(lat)
  contains
    procedure :: init            => tidal_set_init
    procedure :: acceleration    => tide_pressure_accel
    procedure :: get             => get_constituent
  end type

contains

    subroutine tidal_set_init(self, cfg, lat, lon)
      !! Read tidal parameters and initialize derived quantities.
      class(TidalSet),  intent(inout) :: self
      type(ConfigParams), intent(in)  :: cfg
      real(rk), intent(in)            :: lat, lon

      logical :: from_file
      ! Read tidal parameters
      from_file = cfg%get_param_logical('tides.read_from_file', default=.false., strict=.false.)

      if (from_file) then
        call read_from_file(cfg, lat, lon, self)
      else
        call read_from_yaml(cfg, self)
      end if
      ! Create tidal set
      call create_tidal_set(self, lat)
    end subroutine tidal_set_init

    !-------------------------------------------------------------------------
    ! Sum barotropic pressure-gradient accelerations from all tidal constituents.
    !-------------------------------------------------------------------------
    subroutine tide_pressure_accel(self, time, Pxsum, Pysum)
      
      class(TidalSet), intent(in) :: self
      real(rk),             intent(in)  :: time       ! model time since simulation start [s] 
      real(rk),             intent(out) :: Pxsum      ! x and y sea surface slopes for tides in the equation of motion
      real(rk),             intent(out) :: Pysum      ! 
      integer :: i, n
      real(rk) :: sx, sy

      Pxsum = 0.0_rk   ! Here Pxsum and Pysum already sum the accelerations of all constituents, to later call EQN_PRESSURE once
      Pysum = 0.0_rk

      if (.not. allocated(self%C)) return
      n  = size(self%C); if (n == 0) return

      do i = 1, n
        ! Unrotated constituent accelerations
        sx  = - self%C(i)%slopex * cos(self%C(i)%omega * time - self%C(i)%phase)
        sy  = - self%C(i)%slopey * sin(self%C(i)%omega * time - self%C(i)%phase)
        ! Rotate by ellipse orientation (theta) into East/North
          ! Next lines rotate slope forcing to account for tidal ellipse orientation
          !angle = atan2(sy, sx) + self%C(i)%theta
          !Pr = sqrt(sx**2 + sy**2);
        Pxsum = Pxsum + ((self%C(i)%ctheta * sx) - (self%C(i)%stheta * sy)) ! rotate by ellipse orientation to the East
        Pysum = Pysum + ((self%C(i)%stheta * sx) + (self%C(i)%ctheta * sy)) ! rotate by theta to the North
      end do
    end subroutine tide_pressure_accel

    !--------------------------------------------------------------
    ! TidalSet Method to access each constituent values
    !--------------------------------------------------------------
    pure subroutine get_constituent(self, name, c)
      class(TidalSet),  intent(in)  :: self
      character(len=*),    intent(in) :: name
      type(TidalConstituent),   intent(out):: c
      integer :: i
      character(len=:), allocatable :: key

      
      key = to_upper(trim(adjustl(name)))
      c = TidalConstituent()  

      if (.not. allocated(self%c)) return
      do i = 1, size(self%c)
        if (to_upper(trim(self%c(i)%name)) == key) then
          c = self%c(i)
          return
        end if
      end do
    end subroutine get_constituent

    !===========================================================
    !       Internal
    !===========================================================

    subroutine create_tidal_set(S, lat, ok, errmsg)
      !!  - convert smaj/smin from cm s^-1 -> m s^-1
      !!  - convert theta/phase from degrees -> radians
      !!  - assign omega by name
      !!  - compute Coriolis parameter
      real(rk),       intent(in)    :: lat
      type(TidalSet), intent(inout) :: S
      logical,        intent(out), optional :: ok
      character(*),   intent(out), optional :: errmsg

      integer :: i, n
      real(rk) :: w

      if (present(ok))     ok = .true.
      if (present(errmsg)) errmsg = ''

      if (.not. allocated(S%C)) return
      n = size(S%C); if (n == 0) return

      ! Compute coriolis parameter
      S%f0 = 2.0_rk * Omega * sin(deg2rad(lat))

      do i = 1, n
        ! Units: cm/s -> m/s
        S%C(i)%smaj  = 0.01_rk * S%C(i)%smaj
        S%C(i)%smin  = 0.01_rk * S%C(i)%smin

        !Angles: degrees -> radians
        S%C(i)%theta = deg2rad(S%C(i)%theta)
        S%C(i)%phase = deg2rad(S%C(i)%phase)

        ! Assign Omega
        select case (to_upper(S%C(i)%name))
          case ('M2'); w = 1.405278e-4_rk
          case ('S2'); w = 1.454440e-4_rk
          case ('N2'); w = 1.378616e-4_rk
          case ('K1'); w = 7.292952e-5_rk
          case ('O1'); w = 6.759774e-5_rk
          case default
            w = 0.0_rk
        end select
        S%C(i)%omega = w
        if (abs(S%C(i)%smaj) > tiny(1.0_rk)) then
          S%C(i)%polarisn = S%C(i)%smin / S%C(i)%smaj
        else
          S%C(i)%polarisn = 0._rk
        end if
        ! Accelerations in ellipse (major/minor)
        S%C(i)%slopex = (S%C(i)%omega + (S%C(i)%polarisn * S%f0)) * S%C(i)%smaj
        S%C(i)%slopey = (S%f0 + (S%C(i)%polarisn * S%C(i)%omega)) * S%C(i)%smaj
        ! Cos and Sin terms for the rotation matrix ( to rotate by ellipse orientation (theta) into East/North)
        S%C(i)%stheta = sin(S%C(i)%theta)
        S%C(i)%ctheta = cos(S%C(i)%theta)
      end do
      ! Print constituents      
      call print_tides(S)
    end subroutine create_tidal_set  


    subroutine print_tides(tide_set)
      type(TidalSet), intent(in) :: tide_set
      integer  :: i
      print *, '-------------------------------------------------------'
      print *, 'Tidal Constituents loaded: ', size(tide_set%c)
      print *, '-------------------------------------------------------'
      do i = 1, size(tide_set%C)
        write(*,'(A3,2X,"SEMA=",F6.3,2X,"SEMI=",F6.3,2X,"INC=",F5.2, &
              " rad",2X,"PHA=",F5.2," rad")')  &
            trim(tide_set%C(i)%name), tide_set%C(i)%smaj, tide_set%C(i)%smin, &
            tide_set%C(i)%theta, tide_set%C(i)%phase
      end do
      print *, '-------------------------------------------------------'
    end subroutine print_tides


    !---------------- Read tidal parameters from the yaml file ----------------------
    subroutine read_from_yaml(cfg, tide_set)
      implicit none
      type(ConfigParams), intent(in)  :: cfg
      type(TidalSet),  intent(out) :: tide_set

      ! Supported list (extend as needed)
      character(len=*), parameter :: constituents(5) = ['M2','S2','K1','O1','N2']

      integer :: i, nall
      character(len=:), allocatable :: base
      real(rk) :: a, b, ai, ap
      logical  :: fa, fb, found   ! found flags

      nall = size(constituents)
      allocate(tide_set%c(nall))

      do i = 1, nall
        base = 'tides.ellipse_params.' // trim(constituents(i))

        ! Retrieve cnstituents
        a  = cfg%get_param_num(base // '.semi_major', found=fa, nonnegative=.true., finite=.true.)
        b  = cfg%get_param_num(base // '.semi_minor', found=fb, finite=.true.)
        ai = cfg%get_param_num(base // '.inclination', default=0.0_rk, finite=.true.)
        ap = cfg%get_param_num(base // '.phase_ang',   default=0.0_rk, finite=.true.)

        found = fa .and. fb

        if (found) then
          if (abs(b) > a + 1.0e-12_rk) then
            error stop 'Tidal ellipse invalid at '//trim(base)//': |semi_minor| > semi_major.'
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

    subroutine read_from_file(cfg, lat, lon, tide_set)
      type(ConfigParams), intent(in)  :: cfg
      real(rk),           intent(in)  :: lat, lon
      type(TidalSet),     intent(out) :: tide_set

      type(TextTable) :: table

      character(len=*), parameter :: constituents(5) = ['M2','S2','K1','O1','N2']

      character(len=:), allocatable :: fname
      real(rk) :: tol
      real(rk), allocatable :: file_lon(:), file_lat(:), dist(:)
      real(rk), allocatable :: col_values(:)

      integer :: i, k, idx
      integer :: i_lon, i_lat
      integer :: i_sema, i_semi, i_inc, i_pha, nfound

      logical :: ok
      character(len=512) :: errmsg

      fname = cfg%get_param_str('tides.filename', required=.true., trim_value=.true., empty_ok=.false.)
      tol   = cfg%get_param_num('tides.tol_deg', default=0.1_rk, nonnegative=.true., max=10.0_rk, finite=.true.)

      call parse_text_table(fname, table)

      i_lon = table%find_column('lon')
      i_lat = table%find_column('lat')

      if (i_lon == 0) error stop 'read_from_file: lon column not found in tidal file.'
      if (i_lat == 0) error stop 'read_from_file: lat column not found in tidal file.'

      call str_to_real_vec(table%values(:,i_lon), file_lon, ok, errmsg, label='lon')
      if (.not. ok) error stop trim(errmsg)

      call str_to_real_vec(table%values(:,i_lat), file_lat, ok, errmsg, label='lat')
      if (.not. ok) error stop trim(errmsg)

      allocate(dist(table%nrow))

      do i = 1, table%nrow
        dist(i) = simple_distance_deg(lat, lon, file_lat(i), file_lon(i))
      end do

      idx = argmin_abs_vec(dist)

      if (dist(idx) > tol) then
        write(*,'(A,F10.4,A,F10.4,A,F10.4)') &
              'Error reading tidal parameters. Nearest tidal point is too far away: lon=', file_lon(idx), &
              ', lat=', file_lat(idx), ', distance_deg=', dist(idx)
        error stop 'read_from_file: nearest tidal point is outside tides.tol_deg.'
      end if

      allocate(tide_set%C(size(constituents)))

      do k = 1, size(constituents)
        tide_set%C(k)%name  = to_upper(trim(constituents(k)))
        tide_set%C(k)%smaj  = 0.0_rk
        tide_set%C(k)%smin  = 0.0_rk
        tide_set%C(k)%theta = 0.0_rk
        tide_set%C(k)%phase = 0.0_rk

        call get_tidal_column_indices(table, constituents(k), i_sema, i_semi, i_inc, i_pha, nfound)

        if (nfound == 0) then
          write(*,'(A)') 'Warning [read_from_file]: Constituent ' // trim(constituents(k)) // &
                        ' not found in tidal file; set to zero.'
          cycle
        end if

        if (nfound < 4) then
          error stop 'read_from_file: incomplete tidal constituent columns for ' // trim(constituents(k)) // '.'
        end if

        call str_to_real_vec(table%values(:,i_sema), col_values, ok, errmsg, label=trim(table%header(i_sema)))
        if (.not. ok) error stop trim(errmsg)
        tide_set%C(k)%smaj = col_values(idx)

        call str_to_real_vec(table%values(:,i_semi), col_values, ok, errmsg, label=trim(table%header(i_semi)))
        if (.not. ok) error stop trim(errmsg)
        tide_set%C(k)%smin = col_values(idx)

        call str_to_real_vec(table%values(:,i_inc), col_values, ok, errmsg, label=trim(table%header(i_inc)))
        if (.not. ok) error stop trim(errmsg)
        tide_set%C(k)%theta = col_values(idx)

        call str_to_real_vec(table%values(:,i_pha), col_values, ok, errmsg, label=trim(table%header(i_pha)))
        if (.not. ok) error stop trim(errmsg)
        tide_set%C(k)%phase = col_values(idx)

        if (abs(tide_set%C(k)%smin) > tide_set%C(k)%smaj + 1.0e-12_rk) then
          error stop 'read_from_file: invalid tidal ellipse: |semi_minor| > semi_major.'
        end if

      end do

      call table%free()

    end subroutine read_from_file


    subroutine get_tidal_column_indices(table, constituent, i_sema, i_semi, i_inc, i_pha, nfound)
      type(TextTable), intent(in) :: table
      character(*),    intent(in) :: constituent
      integer,         intent(out) :: i_sema, i_semi, i_inc, i_pha
      integer,         intent(out) :: nfound

      character(len=:), allocatable :: cname

      cname = to_upper(trim(constituent))

      i_sema = table%find_column(cname // '_SEMA')
      i_semi = table%find_column(cname // '_SEMI')
      i_inc  = table%find_column(cname // '_INC_DEG')
      i_pha  = table%find_column(cname // '_PHA_DEG')

      nfound = count([i_sema > 0, i_semi > 0, i_inc > 0, i_pha > 0])

    end subroutine get_tidal_column_indices

end module tidal
