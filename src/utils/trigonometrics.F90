module trigonometrics
    use precision_types,  only: rk

    implicit none
    private

    public :: pi, deg2rad
    
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk

contains

  pure function deg2rad(x) result(y)
    real(rk), intent(in) :: x
    real(rk) :: y
    y = x * pi / 180.0_rk
  end function

end module trigonometrics