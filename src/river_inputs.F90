#include "fabm_driver.h"

module pisces_river_inputs
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_river_inputs
      type (type_dependency_id)          :: id_zrivdin
      type (type_state_variable_id)              :: id_tal
      type (type_diagnostic_variable_id) :: id_zrivdin_diag
   contains
       procedure :: initialize
       procedure :: do
    end type

contains
   subroutine initialize(self, configunit)
       class (type_pisces_river_inputs), intent(inout), target :: self
       integer,                  intent(in)            :: configunit


       call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)

       call self%register_dependency(self%id_zrivdin, type_interior_standard_variable(name='input_river', units='molC m-3 s-1'))

       call self%register_diagnostic_variable(self%id_zrivdin_diag, 'zrivdin_diag','-','diagnostic of zrivdin')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_river_inputs), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: zrivdin


      _LOOP_BEGIN_


                _GET_(self%id_zrivdin, zrivdin) ! if ln_trc_cbc(jpno3)=false , zrivdin = 0

                 _ADD_SOURCE_(self%id_tal, - rno3 * zrivdin)

                _SET_DIAGNOSTIC_(self%id_zrivdin_diag, zrivdin) 



      _LOOP_END_

   end subroutine



end module
