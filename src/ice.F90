#include "fabm_driver.h"

module pisces_ice
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_ice
      type (type_surface_dependency_id)          :: id_fmmflx
      type (type_surface_diagnostic_variable_id) :: id_fmmflx_diag, id_zironice
      type (type_state_variable_id)              :: id_fer
      real(rk) :: icefeinput
      logical :: ln_ironice
   contains
       procedure :: initialize
       procedure :: do_surface
    end type

contains
   subroutine initialize(self, configunit)
       class (type_pisces_ice), intent(inout), target :: self
       integer,                  intent(in)            :: configunit

       call self%get_parameter(self%icefeinput, 'icefeinput', 'mol Fe L-1', 'Iron concentration in sea ice', default=15.e-9_rk)
       call self%get_parameter(self%ln_ironice, 'ln_ironice', '-', 'boolean to activate iron inputs from sea ice', default=.false.)
       
       call self%register_dependency(self%id_fmmflx,fmmflx)

       call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')

       call self%register_diagnostic_variable(self%id_fmmflx_diag, 'sfmmflx_diag', 'kg m-2 s-1', 'freshwater budget' )
       call self%register_diagnostic_variable(self%id_zironice, 'Ironice', 'mol Fe m-2 s-1' , 'Iron input/uptake due to sea ice')


   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_ice), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: fer
      real(rk) :: sfmmflx
      real(rk) :: zwflux, zironice


      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_fmmflx, sfmmflx)
         _GET_(self%id_fer, fer)

          zwflux   = sfmmflx / 1000._rk  ! Mokrane: We divide the flux by 1000 to convert from "mol Fe m-3 s-1" to "mol Fe L-1 s-1"
          zironice =  MAX( -0.99_rk * fer, -zwflux * self%icefeinput )

         IF(self%ln_ironice) THEN
           _ADD_SURFACE_FLUX_(self%id_fer, zironice)
         ENDIF

         _SET_SURFACE_DIAGNOSTIC_(self%id_zironice, MAX( -0.99_rk * fer, -zwflux * self%icefeinput * 1.e3))
         _SET_SURFACE_DIAGNOSTIC_(self%id_fmmflx_diag , sfmmflx )

      _SURFACE_LOOP_END_

   end subroutine



end module
