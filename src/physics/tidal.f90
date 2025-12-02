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
  use precision_types,  only: rk
  use physics_params,   only: Omega
  use trigonometrics,   only: pi, deg2rad
  use str_utils,        only: to_upper
  use read_config_yaml, only: ConfigParams
  implicit none
  private
  public :: TidalConstituent, TidalSet 
  public :: create_tidal_set, tide_pressure_slopes
  

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
    procedure :: get => get_constituent
  end type

contains

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


    subroutine tide_pressure_slopes(S, time, Pxsum, Pysum)
      !! Sum barotropic pressure-gradient accelerations from all tidal constituents.
      type(TidalSet),       intent(in)  :: S
      real(rk),             intent(in)  :: time       ! model time since simulation start [s] 
      real(rk),             intent(out) :: Pxsum      ! x and y sea surface slopes for tides in the equation of motion
      real(rk),             intent(out) :: Pysum      ! 
      integer :: i, n
      real(rk) :: sx, sy

      Pxsum = 0.0_rk   ! Here Pxsum and Pysum already sum the accelerations of all constituents, to later call EQN_PRESSURE once
      Pysum = 0.0_rk
      n  = size(S%C); if (n == 0) return

      do i = 1, n
        ! Unrotated constituent accelerations
        sx  = - S%C(i)%slopex * cos(S%C(i)%omega * time - S%C(i)%phase)
        sy  = - S%C(i)%slopey * sin(S%C(i)%omega * time - S%C(i)%phase)
        ! Rotate by ellipse orientation (theta) into East/North
          ! Next lines rotate slope forcing to account for tidal ellipse orientation
          !angle = atan2(sy, sx) + S%C(i)%theta
          !Pr = sqrt(sx**2 + sy**2);
        Pxsum = Pxsum + ((S%C(i)%ctheta * sx) - (S%C(i)%stheta * sy)) ! rotate by ellipse orientation to the East
        Pysum = Pysum + ((S%C(i)%stheta * sx) + (S%C(i)%ctheta * sy)) ! rotate by theta to the North
      end do
    end subroutine  

  !--------------------------------------------------------------
  !----------- TidalSet Method to access each constituent values
  pure subroutine get_constituent(self, name, c)
    class(TidalSet),  intent(in)  :: self
    character(len=*),    intent(in) :: name
    type(TidalConstituent),   intent(out):: c
    integer :: i
    character(len=:), allocatable :: key

    key = to_upper(trim(adjustl(name)))
    c = TidalConstituent()  

    do i = 1, size(self%c)
      if (to_upper(trim(self%c(i)%name)) == key) then
        c = self%c(i)
        return
      end if
    end do
  end subroutine get_constituent

  subroutine print_tides(tide_set)
    type(TidalSet), intent(in) :: tide_set
    integer  :: i
    print *, '-------------------------------------------------------'
    print *, 'Tidal Constitituents loaded: ', size(tide_set%c)
    print *, '-------------------------------------------------------'
    do i = 1, size(tide_set%C)
      write(*,'(A3,2X,"SEMA=",F6.3,2X,"SEMI=",F6.3,2X,"INC=",F5.2, &
             " rad",2X,"PHA=",F5.2," rad")')  &
          trim(tide_set%C(i)%name), tide_set%C(i)%smaj, tide_set%C(i)%smin, &
          tide_set%C(i)%theta, tide_set%C(i)%phase
    end do
    print *, '-------------------------------------------------------'
  end subroutine print_tides

end module tidal
