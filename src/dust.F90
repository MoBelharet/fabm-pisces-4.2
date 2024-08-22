#include "fabm_driver.h"

module pisces_dust
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_dust
      type (type_surface_dependency_id)          :: id_dustdep, id_fmmflx
      type (type_dependency_id)                  :: id_gdept_n
      type (type_diagnostic_variable_id)         :: id_zdust
      type (type_surface_diagnostic_variable_id) :: id_zirondep, id_pdust
      type (type_state_variable_id)              :: id_sil, id_po4, id_fer
      logical :: ln_ironice
      real(rk) ::  mfrac, wdust, icefeinput
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_dust), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_implemented_routines((/source_do, source_do_surface/))

      call self%get_parameter(self%ln_ironice, 'ln_ironice','', 'boolean for Fe input from sea ice', default=.false.)
      call self%get_parameter(self%icefeinput, 'icefeinput', 'mol Fe L-1', 'Iron concentration in sea ice', default=15.e-9_rk)
      call self%get_parameter(self%mfrac, 'mfrac', '1', 'Fe mineral fraction of dust', default=0.035_rk)
      call self%get_parameter(self%wdust, 'wdust', 'm d-1', 'sinking speed of dust', default=2._rk)

      call self%register_dependency(self%id_dustdep,dustdep) ! Mokrane dustdep is declared in shared.F90 file
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_fmmflx,fmmflx)

      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_sil, 'sil', 'mol Si L-1', 'silicate')

      call self%register_diagnostic_variable(self%id_zirondep, 'zirondep', 'mol m-2 s-1', 'iron deposition')
      call self%register_diagnostic_variable(self%id_pdust, 'pdust', 'g m-3', 'concentration at the surface')
      call self%register_diagnostic_variable(self%id_zdust, 'zdust', 'g m-3', 'concentration')
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_pisces_dust), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: fer
      real(rk) :: sfmmflx
      real(rk) :: zwflux, zironice


      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_fmmflx, sfmmflx)
         _GET_(self%id_fer, fer)

         zwflux   = sfmmflx / 1000._rk  ! Mokrane: We divide the flux by 1000 to convert from "mol Fe m-3 s-1" to "mol Fe L-1 s-1"
         zironice =  MAX( -0.99_rk * fer, -zwflux * self%icefeinput )

        IF(self%ln_ironice) _ADD_SURFACE_FLUX_(self%id_fer, zironice)
         
      _SURFACE_LOOP_END_
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_dust), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: dust, gdept_n, zdust, zirondep, zpdep, zdustdep, zwdust

      _LOOP_BEGIN_
         _GET_SURFACE_(self%id_dustdep, dust)
         _GET_(self%id_gdept_n, gdept_n)
         
         !************ Mokrane **********
         zwdust = 0.03_rk / ( self%wdust / rday ) / ( 250._rk * rday )
        zdustdep = dust * zwdust  * EXP( -gdept_n /( 250._rk * self%wdust ) ) !* rfact
        _ADD_SOURCE_(self%id_fer, + zdustdep * self%mfrac / mMass_Fe)

        _ADD_SOURCE_(self%id_po4 ,  + zdustdep * 1.e-3 / mMass_P)

        _ADD_SOURCE_(self%id_sil,   + zdustdep * 0.269_rk / mMass_Si)

        
        !************************************


         zdust = dust / ( self%wdust / rday )
         _SET_DIAGNOSTIC_(self%id_zdust, zdust)

      _LOOP_END_
   end subroutine

end module
