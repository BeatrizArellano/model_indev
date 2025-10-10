module precision_types
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, real128
  implicit none
  private                          

  ! ---- Integers ----
  integer, parameter, public :: ik = int32    ! default integer kind
  integer, parameter, public :: lk = int64    ! large counters/indices

  ! ---- Reals ----
  integer, parameter, public :: rk = real64   ! model real kind (default)
  integer, parameter, public :: qk = real128  ! high-precision kind (quad)

end module precision_types
