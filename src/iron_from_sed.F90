#include "fabm_driver.h"

module pisces_iron_from_sed
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_iron_from_sed
      type (type_state_variable_id)              :: id_fer
      type (type_dependency_id)                  :: id_zcmask_cal, id_zcmask_ext, id_gdepw_n, id_e3t_0
      type (type_diagnostic_variable_id)         :: id_ironsed
      real(rk) :: sedfeinput
   contains
      procedure :: initialize
      procedure :: do
   end type
contains
   subroutine initialize(self, configunit)
        class (type_pisces_iron_from_sed), intent(inout), target :: self
        integer,                      intent(in)            :: configunit

        call self%get_parameter(self%sedfeinput, 'sedfeinput', 'mol Fe L-1 d-1 m', 'iron flux from the sediments', default=2.e-9_rk)

        call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')

        call self%register_dependency(self%id_zcmask_cal, type_interior_standard_variable(name='coastal_island_mask', units='1'))
        call self%register_dependency(self%id_e3t_0, type_interior_standard_variable(name='ref_cell_thickness' , units='m') )
        call self%register_dependency(self%id_zcmask_ext, bathy_etop5)

        call self%register_dependency(self%id_gdepw_n, standard_variables%depth)
        
        call self%register_diagnostic_variable(self%id_ironsed, 'ironsed', 'mol Fe m-2 s-1',  'iron inputs')

   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
        class (type_pisces_iron_from_sed), intent(in) :: self

          _DECLARE_ARGUMENTS_DO_
             real(rk) :: zcmask_cal, zcmask_ext, gdepw_n
             real(rk) :: zcmask, zexpide, zdenitide, ironsed, e3t_0

          _LOOP_BEGIN_
                _GET_(self%id_zcmask_cal, zcmask_cal)
                _GET_(self%id_zcmask_ext, zcmask_ext)
                _GET_(self%id_gdepw_n, gdepw_n)
                _GET_(self%id_e3t_0, e3t_0)

                zcmask = zcmask_cal + zcmask_ext

                zexpide   = MIN( 8._rk,( gdepw_n / 500. )**(-1.5) )                           ! Eq 85a - taken from p4zsbc.F90 but replaced cell center depth with bottom depth
                zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2  ! Eq 85b
                ironsed = zcmask * MIN( 1._rk, EXP( zdenitide ) / 0.5 )                       ! Eq 85c
                ironsed = self%sedfeinput * ironsed * r1_rday / e3t_0

                _SET_DIAGNOSTIC_(self%id_ironsed, ironsed)
                _ADD_SOURCE_(self%id_fer, ironsed)



          _LOOP_END_

  end subroutine


end module







