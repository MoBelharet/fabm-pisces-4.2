#include "fabm_driver.h"

module pisces_size_phyto
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_size_phyto
      type (type_diagnostic_variable_id) :: id_sizep_prev
      type (type_dependency_id)          :: id_sizep
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_size_phyto), intent(inout), target :: self
      integer,                   intent(in)            :: configunit


      call self%register_diagnostic_variable(self%id_sizep_prev, 'sizep_prev','-','Mean relative size')
      call self%register_dependency(self%id_sizep, 'sizep','-','Mean relative size')


   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_size_phyto), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sizep, sizep_prev

      _LOOP_BEGIN_
         _GET_(self%id_sizep, sizep)

         sizep_prev = sizep

         _SET_DIAGNOSTIC_(self%id_sizep_prev, sizep_prev)
      _LOOP_END_
   end subroutine

end module
