module freshwater_fluxes
   use precision_types, only: rk
   implicit none
   private

   public :: apply_surface_freshwater

contains


   !!  Subroutine to compute surface freshwater fluxes and salinity changes
   !-------------------------------------------------------------------
   subroutine apply_surface_freshwater(salt, dz, precip, evap, runoff)
      real(rk), intent(inout) :: salt(:)        ! Salinity (vertical profile)
      real(rk), intent(in)    :: dz(:)          ! Layer thicknesses
      real(rk), intent(in)    :: precip, evap   ! In m/s 
      real(rk), intent(in), optional :: runoff

      ! Evaporation in ERA5 is negative for evaporation and positive for condensation

      !write(*,'(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
      !      'precip = ', precip, 'evap = ',   evap, 'runoff = ', runoff
      return
   end subroutine apply_surface_freshwater

end module freshwater_fluxes
