#include "fabm_driver.h"

module pisces_island_coastal_masks
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_island_coastal_masks
      type (type_dependency_id)                  :: id_zcmask_cal, id_zcmask_ext
      type (type_diagnostic_variable_id)         :: id_zcmask
   contains
      procedure :: initialize
      procedure :: do
   end type
contains
   subroutine initialize(self, configunit)
        class (type_pisces_island_coastal_masks), intent(inout), target :: self
        integer,                      intent(in)            :: configunit

        call self%register_dependency(self%id_zcmask_cal, type_interior_standard_variable(name='coastal_island_mask', units='1'))
        call self%register_dependency(self%id_zcmask_ext, bathy_etop5)

        call self%register_diagnostic_variable(self%id_zcmask, 'zcmask', '1', 'Fractional area for the bottom bathymetry')

   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
        class (type_pisces_island_coastal_masks), intent(in) :: self

          _DECLARE_ARGUMENTS_DO_
             real(rk) :: zcmask_cal, zcmask_ext

          _LOOP_BEGIN_
                _GET_(self%id_zcmask_cal, zcmask_cal)
                _GET_(self%id_zcmask_ext, zcmask_ext)

                _SET_DIAGNOSTIC_(self%id_zcmask, zcmask_cal + zcmask_ext)



          _LOOP_END_

  end subroutine


end module







