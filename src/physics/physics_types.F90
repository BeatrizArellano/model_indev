module physics_types
    use precision_types,     only: rk, lk
    use grids,               only: VerticalGrid
    use tidal,               only: TidalSet
    use tridiagonal,         only: TridiagCoeff
    use physics_params,      only: PhysicsParams, z0s_min

  implicit none

  public :: PhysicsState, PhysicsEnv, Dirichlet, Neumann

  !======================
  ! Internal state
  !======================

  type :: PhysicsState
    integer              :: N                     ! Number of layers in the water column
    ! --- Prognostic variables at layers' centres (1..N) ---
    real(rk), allocatable :: temp(:), sal(:), rho(:)  ! temperature, salinity and density
    real(rk), allocatable :: velx(:), vely(:)         ! Velocity components

    ! --- Turbulence and mixing on layer interfaces (0..N) ---
    real(rk), allocatable :: Kz(:)                    ! scalar diffusivity at interfaces [m2/s]
    real(rk), allocatable :: Nz(:)                    ! momentum viscosity at faces [m2/s], 0..N
    real(rk), allocatable :: tke(:)                   ! TKE at interfaces
    real(rk), allocatable :: eps(:)                   ! epsilon
    real(rk), allocatable :: Lscale(:)                ! Length scale
    real(rk), allocatable :: cmue1(:)                 ! Stability function 1 (used outside to pass to dissipation)

    ! --- Surface and bottom ---
    real(rk) :: tau_x = 0.0_rk, tau_y = 0.0_rk        ! wind stress components [N m-2]
    real(rk) :: u_taus = 0.0_rk, u_taub = 0.0_rk      ! friction velocities at surface and bottom [m s-1]
    real(rk) :: z0s    = z0s_min, z0b  = 0.01_rk      ! surface and bottom roughness lengths [m]

    ! --- Heat fluxes  ---
    !real(rk) :: Q_sw_net = 0.0_rk, Q_lw_net = 0.0_rk, Q_lat = 0.0_rk, Q_sens = 0.0_rk, Q_net = 0.0_rk
  end type PhysicsState

  ! An envelope for the Physics Environment
  type :: PhysicsEnv   
    type(VerticalGrid)   :: grid       ! per-column grid (copy or pointer to external)
    type(PhysicsState)   :: PS         ! prognostic + turbulence arrays
    type(PhysicsParams)  :: params     ! per-column physics params (or shared elsewhere)
    type(TidalSet)       :: Tides      ! site-specific tides
    type(TridiagCoeff)   :: trid        ! workspace for implicit solves
    logical              :: is_init = .false. ! Flag to indicate that the physics is initialised for this water column
  end type PhysicsEnv

  ! Boundary types
  integer, parameter :: Dirichlet = 1
  integer, parameter :: Neumann   = 2
end module physics_types
