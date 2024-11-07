#include "fabm_driver.h"

module pisces_nutrient_deposition
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_nutrient_deposition
      type (type_surface_dependency_id)          :: id_zndep_no3, id_zndep_nh4
      type (type_state_variable_id)              :: id_tal
      type (type_surface_diagnostic_variable_id) :: id_zndep_no3_diag, id_zndep_nh4_diag
      type (type_dependency_id)                  :: id_e3t_n
   contains
       procedure :: initialize
       procedure :: do_surface
    end type

contains
   subroutine initialize(self, configunit)
       class (type_pisces_nutrient_deposition), intent(inout), target :: self
       integer,                  intent(in)            :: configunit


       call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)

       call self%register_dependency(self%id_zndep_no3, type_horizontal_standard_variable(name='NO3_deposition_flux', units='1'))
       call self%register_dependency(self%id_zndep_nh4, type_horizontal_standard_variable(name='NH4_deposition_flux', units='1'))

       call self%register_diagnostic_variable(self%id_zndep_no3_diag, 'zndep_no3_diag','-','diagnostic of zndep_no3')
       call self%register_diagnostic_variable(self%id_zndep_nh4_diag, 'zndep_nh4_diag','-','diagnostic of zndep_nh4')

       call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)

   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_nutrient_deposition), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: tal
      real(rk) :: zndep_no3, zndep_nh4, e3t_n

       
      _SURFACE_LOOP_BEGIN_


      _GET_SURFACE_(self%id_zndep_no3, zndep_no3) ! if ln_trc_sbc(jpno3)=false , zndep_no3 = 0
      _GET_SURFACE_(self%id_zndep_nh4, zndep_nh4) ! if ln_trc_sbc(jpnh4)=false , zndep_nh4 = 0
      !_GET_(self%id_tal, tal) !
      _GET_(self%id_e3t_n, e3t_n)

        _ADD_SURFACE_FLUX_(self%id_tal, - rno3 * zndep_no3)

        _ADD_SURFACE_FLUX_(self%id_tal, + rno3 * zndep_nh4)

        _SET_SURFACE_DIAGNOSTIC_(self%id_zndep_no3_diag, zndep_no3/ e3t_n) 
        _SET_SURFACE_DIAGNOSTIC_(self%id_zndep_nh4_diag, zndep_nh4 / e3t_n ) 


      _SURFACE_LOOP_END_

   end subroutine



end module
