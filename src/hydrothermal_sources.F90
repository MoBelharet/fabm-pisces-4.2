#include "fabm_driver.h"

module pisces_hydrothermal_sources
   use fabm_types
   use pisces_shared

   implicit none
   
   private

   type, extends(type_base_model), public :: type_pisces_hydrothermal_sources
      type (type_surface_dependency_id)          :: id_cell_area
      type (type_dependency_id)                  :: id_e3t_n, id_hydrofe
      type (type_state_variable_id)              :: id_fer
      type (type_diagnostic_variable_id)         :: id_hydrofe_diag
      real(rk) :: hratio
!      logical  :: ln_hydrofe
   contains
      procedure :: initialize
      procedure :: do
   end type
contains

   subroutine initialize(self, configunit)
        class (type_pisces_hydrothermal_sources), intent(inout), target :: self
        integer,                      intent(in)            :: configunit

        call self%register_implemented_routines((/source_do/))

        call self%get_parameter(self%hratio, 'hratio', '1', 'Fe to 3He ratio assumed for vent iron supply', default=1.e+7_rk)
        !      call self%get_parameter(self%ln_hydrofe, 'ln_hydrofe','', 'boolean for Fe input from sea hydrothermal sources', default=.false.)
        call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
        call self%register_dependency(self%id_hydrofe, iron_hydrothermal_source)
        !call self%register_dependency(self%id_cell_area, cell_area)
        call self%register_dependency(self%id_cell_area, standard_variables%cell_area)
        call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)

        call self%register_diagnostic_variable(self%id_hydrofe_diag, 'hydrofe_diag','mol Fe L-1 s-1', 'diagnostic of hydrofe')

   end subroutine


   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_hydrothermal_sources), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: carea, e3t_n, sf_hydrofe, hydrofe
      real(rk), parameter :: ryyss    = nyear_len * rday    ! number of seconds per year, Jorn: nyear_len should account for leap years, calendars, etc.

      _LOOP_BEGIN_

          _GET_SURFACE_(self%id_cell_area, carea)
          _GET_(self%id_e3t_n, e3t_n)
          _GET_(self%id_hydrofe, sf_hydrofe)

          hydrofe = ( MAX( rtrn , sf_hydrofe ) * self%hratio ) &
           &       / ( carea * e3t_n * ryyss + rtrn ) / 1000._rk

         _SET_DIAGNOSTIC_(self%id_hydrofe_diag, hydrofe * 1.e+3)

         _ADD_SOURCE_(self%id_fer, + hydrofe)


      _LOOP_END_

   end subroutine


end module
















