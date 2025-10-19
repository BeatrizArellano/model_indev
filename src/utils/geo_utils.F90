! Useful functions for gridpoints and geographical distance
module geo_utils
  use precision_types, only: rk    ! Importing real64 and int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  public :: is_lat_valid, is_lon_valid, is_depth_valid
  public :: convert_to_lon180, convert_to_lon360
  public :: simple_distance_deg

  type, public :: LocationInfo
    character(len=:), allocatable :: name
    real(rk) :: lat
    real(rk) :: lon
    real(rk) :: depth                       ! depth [m]
  end type

contains



    ! Verifies if latitude is finite and within [-90, 90]
    pure logical function is_lat_valid(lat) result(ok)
        real(rk), intent(in) :: lat
        ok = ieee_is_finite(lat) .and. (lat >= -90._rk) .and. (lat <= 90._rk)
    end function is_lat_valid

    ! Verifies if longitude is finite and within [-360,360]
    pure logical function is_lon_valid(lon) result(ok)
        real(rk), intent(in) :: lon
        ok = ieee_is_finite(lon) .and. (lon >= -360._rk) .and. (lon <= 360._rk)
    end function is_lon_valid


    ! Convert any longitude to [-180, 180)
    pure real(rk) function convert_to_lon180(lon) result(lon180)
        real(rk), intent(in) :: lon
        lon180 = modulo(lon + 180._rk, 360._rk) - 180._rk
        if (lon180 == -180._rk) lon180 = 180._rk  ! Convert -180 to +180
    end function convert_to_lon180

    pure real(rk) function convert_to_lon360(lon) result(lon360)
        real(rk), intent(in) :: lon
        lon360 = modulo(modulo(lon,360._rk) + 360._rk, 360._rk)
    end function convert_to_lon360


    ! Verifies that depth is finite and >0
    pure logical function is_depth_valid(depth) result(ok)
        real(rk), intent(in) :: depth
        ok = ieee_is_finite(depth) .and. (depth > 0._rk)
    end function is_depth_valid

    ! Calculates a simple distance in degrees betw two grid-points
    pure function simple_distance_deg(lat1, lon1, lat2, lon2) result(ddeg)
      real(rk), intent(in) :: lat1, lon1, lat2, lon2   ! degrees
      real(rk) :: ddeg, dlat, dlon, cosphi
      real(rk), parameter :: DEG2RAD = acos(-1.0_rk)/180.0_rk

      dlat  = lat2 - lat1
      dlon  = lon2 - lon1
      cosphi = cos( 0.5_rk*(lat1+lat2)*DEG2RAD )
      ddeg = sqrt( dlat*dlat + (dlon*cosphi)*(dlon*cosphi) )
   end function simple_distance_deg

end module geo_utils