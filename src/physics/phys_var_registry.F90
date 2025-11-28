module phys_var_registry
    use read_config_yaml,  only: ConfigParams
    use physics_types,     only: PhysicsEnv
    use variable_registry, only: VarMetadata, register_variable

  implicit none

  public :: register_physics_variables

contains

    subroutine register_physics_variables(cfg_params, PE)
      type(ConfigParams), intent(in)    :: cfg_params
      type(PhysicsEnv),   intent(inout) :: PE

      ! Temperature [°C] 
      call register_variable(PE%phys_vars, name='temperature', standard_name='sea_water_temperature', &
                             long_name='Sea water temperature (in-situ).',                                       &
                             units='degC', vert_coord='centre', n_space_dims=1, data_1d=PE%PS%temp)

      ! Salinity 
      call register_variable(PE%phys_vars, name='salinity', standard_name='sea_water_practical_salinity', &
                             long_name='Sea water practical salinity.',                             &
                             units='PSU', vert_coord='centre', n_space_dims=1, data_1d=PE%PS%sal)

      ! Density [kg m-3] 
      call register_variable(PE%phys_vars, name='density', standard_name='sea_water_density', &
                             long_name='Sea water density (in-situ).',                        &
                             units='kg m-3', vert_coord='centre', n_space_dims=1, data_1d=PE%PS%rho)

      ! Velocity components 
      call register_variable(PE%phys_vars, name='u', standard_name='eastward_sea_water_velocity', &
                             long_name='Eastward seawater velocity from barotropic tidal pressure-gradient accelerations.', &
                             units='m s-1', vert_coord='centre', n_space_dims=1, data_1d=PE%PS%velx)

      call register_variable(PE%phys_vars, name='v', standard_name='northward_sea_water_velocity', &
                             long_name='Northward seawater velocity from barotropic tidal pressure-gradient accelerations.', &
                             units='m s-1', vert_coord='centre', n_space_dims=1, data_1d=PE%PS%vely)

      !-----------------------------------
      ! Variables at layer interfaces
      !-----------------------------------

      ! Vertical eddy diffusivity [m2 s-1]
      call register_variable(PE%phys_vars, name='Kz', standard_name='vertical_eddy_diffusivity_of_tracers', &
                             long_name='Vertical eddy diffusivity.', &
                             units='m2 s-1', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%Kz)

      ! Momentum viscosity [m2 s-1]
      call register_variable(PE%phys_vars, name='Nz', standard_name='ocean_vertical_momentum_diffusivity', &
                             long_name='Vertical eddy viscosity for momentum.',                                   &
                             units='m2 s-1', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%Nz)

      ! TKE [m2 s-2] – CF: specific_turbulent_kinetic_energy_of_sea_water
      call register_variable(PE%phys_vars, name='tke', standard_name='specific_turbulent_kinetic_energy_of_sea_water', &
                             long_name='Turbulent kinetic energy.',                                          &
                             units='m2 s-2', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%tke)

      ! Dissipation rate [m2 s-3] 
      call register_variable(PE%phys_vars, name='eps', standard_name='specific_turbulent_kinetic_energy_dissipation_in_sea_water', &
                             long_name='Turbulent kinetic energy dissipation rate.',                   &
                             units='m2 s-3', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%eps)

      call register_variable(PE%phys_vars, name='NN', &
                             standard_name='square_of_buoyancy_frequency_in_sea_water', &
                             long_name='Squared buoyancy frequency.', &
                             units='s-2', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%NN)

      call register_variable(PE%phys_vars, name='SS', &
                             long_name='Squared shear frequency.', &
                             units='s-2', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%SS)

       call register_variable(PE%phys_vars, name='Ri', standard_name='richardson_number_in_sea_water', &
                              long_name='Gradient Richardson number', &
                              units='', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%Ri)

      ! Turbulent length scale [m] 
      call register_variable(PE%phys_vars, name='Lscale', long_name='Turbulent mixing length scale', &                            
                             units='m', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%Lscale)

      ! Stability function cmue1 (units?)
      call register_variable(PE%phys_vars, name='cmue1', long_name='Stability function for momentum diffusivity in the Canuto turbulence scheme.', &
                             units='', vert_coord='interface', n_space_dims=1, data_1d=PE%PS%cmue1)



      !-----------------------------------
      ! Surface
      !-----------------------------------

      ! Wind stress components at the surface 
      call register_variable(PE%phys_vars, name='tau_x', standard_name='surface_downward_eastward_stress', &                             
                             long_name='Eastward component of surface wind stress', &
                             units='N m-2', vert_coord='surface', n_space_dims=0, data_0d=PE%PS%tau_x)


      call register_variable(PE%phys_vars, name='tau_y', standard_name='surface_downward_northward_stress', &
                             long_name='Northward component of surface wind stress', &
                             units='N m-2', vert_coord='surface', n_space_dims=0, data_0d=PE%PS%tau_y)

      ! Surface friction velocity 
      call register_variable(PE%phys_vars, name='u_taus',                                           &
                             long_name='Water-side surface friction velocity.', &
                             units='m s-1', vert_coord='surface', n_space_dims=0, data_0d=PE%PS%u_taus)

      ! Surface roughness length [m]
      call register_variable(PE%phys_vars, name='z0s', standard_name='surface_roughness_length', &
                             long_name='Water-side surface roughness length.', &
                             units='m', vert_coord='surface', n_space_dims=0, data_0d=PE%PS%z0s)

      !-----------------------------------
      ! Bottom
      !-----------------------------------

      ! Bottom friction velocity [m s-1]
      call register_variable(PE%phys_vars, name='u_taub', standard_name='bottom_friction_velocity', &
                             long_name='Bottom friction velocity.', &
                             units='m s-1', vert_coord='bottom', n_space_dims=0, data_0d=PE%PS%u_taub)

      ! Bottom roughness length [m]
      call register_variable(PE%phys_vars, name='z0b', standard_name='bottom_roughness_length', &
                             long_name='Bottom roughness length.',                  &
                             units='m', vert_coord='bottom', n_space_dims=0, data_0d=PE%PS%z0b)

      ! Bottom stress magnitude [N m-2]
      call register_variable(PE%phys_vars, name='stressb',                           &
                             long_name='Bottom stress magnitude (Pa).',              &
                             units='N m-2', vert_coord='bottom', n_space_dims=0, data_0d=PE%PS%stressb)


       call define_output_variables(cfg_params,PE%phys_vars)
    end subroutine register_physics_variables

    subroutine define_output_variables(cfg_params,VMD)
      type(ConfigParams), intent(in)    :: cfg_params
      type(VarMetadata),   intent(inout) :: VMD(:)

      logical :: output_all, output_var, default
      integer :: i
      character(len=:), allocatable :: key

      output_all = cfg_params%get_param_logical('physics.output.all', default=.false.)

      if (output_all) then
        ! Enable everything
        do i = 1, size(VMD)
            VMD(i)%output = .true.
        end do
        return
      end if

      ! Otherwise, look at which variables to output
      do i = 1, size(VMD)
        ! Build YAML key
        key = 'physics.output.' // trim(VMD(i)%name)
        ! Get values. Default is false, except for temperature
        if (trim(VMD(i)%name) == 'temperature') then
          default = .true.
        else
          default = .false.
        end if
        output_var = cfg_params%get_param_logical(key, default=default)
        VMD(i)%output = output_var
      end do
    end subroutine define_output_variables
end module phys_var_registry
