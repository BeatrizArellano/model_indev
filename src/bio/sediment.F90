module sediment    
    use bio_params,        only: BioParams
    use bio_types,         only: BioEnv
    use fabm,              only: type_fabm_model, type_fabm_interior_variable_id, &
                                   type_fabm_horizontal_variable_id, type_fabm_scalar_variable_id
    use grids,             only: VerticalGrid
    use precision_types,   only: rk, lk
    use read_config_yaml,  only: ConfigParams
    
    use tridiagonal,       only: TridiagCoeff
    use variable_registry, only: VarMetadata

implicit none

contains



    ! Initialise sediments
    subroutine init_sediments(cfg, SedE)
        type(ConfigParams),       intent(in) :: cfg
        type(BioEnv),          intent(inout) :: SedE

    end subroutine init_sediments


    ! Helpers

    


end module sediment