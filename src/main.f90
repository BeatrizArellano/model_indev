program main
  use iso_fortran_env, only: real64
  use read_config_yaml
  use physics_driver
  implicit none

  type(ConfigParams) :: cfg_user
  real(real64) :: dt, t
  integer :: nsteps, n

  call cfg_user%init()
  call cfg_user%load_yaml_content('physics.yaml','physics')
  call cfg_user%load_yaml_content('main.yaml','main')

  call physics_init(cfg_user)

  call physics_end()

end program main




