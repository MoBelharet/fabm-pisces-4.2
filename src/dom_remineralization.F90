#include "fabm_driver.h"

module pisces_dom_remineralization
   use fabm_types
   use pisces_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_pisces_dom_remineralization
      type (type_state_variable_id) :: id_no3, id_nh4, id_po4, id_fer, id_doc, id_oxy, id_dic, id_tal, id_sfe, id_bfe
      type (type_dependency_id) :: id_zoo, id_mes, id_gdept_n, id_tem
      type (type_surface_dependency_id) :: id_hmld, id_heup, id_heup_01
      type (type_dependency_id) :: id_zdepbac, id_nitrfac
      type (type_diagnostic_variable_id) :: id_remin, id_denit
      type (type_diagnostic_variable_id) :: id_febact, id_blim

      real(rk) :: xremikc, xkdoc, concno3, concnh4, concfe
      real(rk) :: feratb, xkferb, mumax0
   contains
      procedure :: initialize
      procedure :: do
   end type

   type, extends(type_base_model), public :: type_pisces_bacteria
      type (type_dependency_id) :: id_gdept_n, id_zoo, id_mes
      type (type_surface_dependency_id) :: id_hmld, id_heup_01
      type (type_diagnostic_variable_id) :: id_bact
   contains
      procedure :: initialize => bacteria_initialize
      procedure :: do_column  => bacteria_do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_dom_remineralization), intent(inout), target :: self
      integer,                                  intent(in)            :: configunit

      class (type_pisces_bacteria), pointer :: pbacteria

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(self%xremikc, 'xremikc', 'd-1', 'remineralization rate of doc', default=0.4_rk)
      call self%get_parameter(self%xkdoc, 'xkdoc', 'mol C L-1', 'DOC half-saturation constant', default=417.E-6_rk)
      call self%get_parameter(self%concno3, 'concno3', 'mol C L-1', 'nitrate half-saturation constant', default=3.E-7_rk)   ! ???? 0.03 umol N/L in paper
      call self%get_parameter(self%concnh4, 'concnh4', 'mol C L-1', 'ammonium/phosphate half-saturation constant', default=3.E-7_rk)! ???? 0.003 umol N or P/L in paper
      call self%get_parameter(self%concfe, 'concfe', 'mol Fe L-1', 'iron half-saturation constant', default=3.E-11_rk)

      call self%get_parameter(self%mumax0, 'mumax0', 'd-1', 'maximum iron uptake rate of bacteria at 0 degrees Celsius', default=0.6_rk)
      call self%get_parameter(self%feratb, 'feratb', 'mol Fe (mol C)-1', 'Fe/C quota in bacteria', default=60.E-6_rk)
      call self%get_parameter(self%xkferb, 'xkferb', 'mol Fe L-1', 'half-saturation constant for bacteria Fe/C', default=4.E-10_rk)

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)

      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup, 'heup', 'm', 'euphotic depth')
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'depth where daily mean PAR equals 0.5 W m-2')
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      !call self%register_dependency(self%id_hmld, standard_variables%mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_dependency(self%id_mes, 'mes', 'mol C L-1', 'mesozooplankton')
      call self%register_dependency(self%id_tem, standard_variables%temperature)

      allocate(pbacteria)
      call self%add_child(pbacteria, 'bacteria')
      call pbacteria%request_coupling(pbacteria%id_heup_01, '../heup_01')
      call pbacteria%request_coupling(pbacteria%id_zoo, '../zoo')
      call pbacteria%request_coupling(pbacteria%id_mes, '../mes')

      call self%register_dependency(self%id_zdepbac, 'zdepbac', 'mol C L-1', 'bacterial biomass proxy')
      call self%request_coupling(self%id_zdepbac, './bacteria/bact')
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')

      call self%register_diagnostic_variable(self%id_remin, 'remin', 'mol N L-1 s-1', 'remineralization')
      call self%register_diagnostic_variable(self%id_denit, 'denit', 'mol N L-1 s-1', 'denitrification')
      call self%register_diagnostic_variable(self%id_febact, 'febact', 'mol Fe L-1 s-1', 'bacterial uptake of Fe')
      call self%register_diagnostic_variable(self%id_blim, 'blim', '-', 'bacterial limitation of remineralization')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_dom_remineralization), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: no3, nh4, po4, fer, doc, oxy
      real(rk) :: hmld, heup, gdept_n, tem, heup_01
      real(rk) :: zdenom, xnanono3, xnanonh4
      real(rk) :: zlim1, zlim2, zlim3, zlim4, xlimbacl, xlimbac
      real(rk) :: zdepbac, nitrfac, nitrfac2
      real(rk) :: zremik, zolimic, zolimi, zammonic, denitr, zoxyremc
      real(rk) :: tgfunc, zdep, zdepmin, zdepprod, zdepeff, zbactfer, blim, zremikc, zlimnh4, zlimno3, znutlimtot !,zoo, mes, ztempbac

      real(rk), parameter :: xstep = r1_rday

      _LOOP_BEGIN_
         _GET_(self%id_no3, no3)
         _GET_(self%id_nh4, nh4)
         _GET_(self%id_po4, po4)
         _GET_(self%id_fer, fer)
         _GET_(self%id_doc, doc)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_zdepbac, zdepbac)
         _GET_(self%id_nitrfac, nitrfac)
         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_SURFACE_(self%id_heup, heup)
         _GET_SURFACE_(self%id_heup_01, heup_01)
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_tem, tem)                  ! temperature (degrees Celsius)

         ! Jorn: from p4zlim.F90
         !zdenom = 1._rk /  ( self%concno3 * self%concnh4 + self%concnh4 * no3 + self%concno3 * nh4 )  ! Jorn: denominator in Eq 34g,h
         !xnanono3 = no3 * self%concnh4 * zdenom                              ! Jorn: Eq 34h
         !xnanonh4 = nh4 * self%concno3 * zdenom                    ! Jorn: Eq 34g
!********************/ Mokrane / ************************
         
         zlimnh4 = nh4 / ( self%concno3 + nh4 )
         zlimno3 = no3 / ( self%concno3 + no3 )
         znutlimtot = ( nh4 + no3 ) / ( self%concno3 + nh4 + no3 )
         xnanonh4 = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
         xnanono3 = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )

!********************************************************

         !
         zlim1    = xnanono3 + xnanonh4                  ! Jorn: Eq 34f
         zlim2    = po4 / ( po4 + self%concnh4 )         ! Jorn: Eq 34e
         zlim3    = fer / ( self%concfe + fer )          ! Jorn: Eq 34d
         zlim4    = doc / ( self%xkdoc   + doc )         ! Jorn: Eq 34b
         xlimbacl = MIN( zlim1, zlim2, zlim3 )           ! Jorn: Eq 34c, nutrient limitation, but uses total nitrogen limitation rather than ammonium-only limitation as in paper (but that likely is a typo)
         xlimbac  = MIN( zlim1, zlim2, zlim3 ) * zlim4   ! Jorn: Eq 34a, nutrient and substrate limitation

         ! denitrification factor computed from NO3 levels
         nitrfac2 = MAX( 0.e0_rk,       ( 1.E-6_rk - no3 )  &
            &                                / ( 1.E-6_rk + no3 ) )
         nitrfac2 = MIN( 1., nitrfac2 )

         ! Jorn: from p4zrem.F90
         ! DOC ammonification. Depends on depth, phytoplankton biomass
         ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 

         zremik = xstep / 1.e-6 * xlimbac * zdepbac
         zremik = MAX( zremik, 2.74e-4 * xstep / self%xremikc )
         zremikc = self%xremikc * zremik
         ! Ammonification in oxic waters with oxygen consumption
         ! -----------------------------------------------------

         zolimic = zremikc * ( 1._rk- nitrfac ) * doc 
         zolimi =  MAX(0._rk, MIN( ( oxy - rtrn ) / o2ut, zolimic ) )

         ! Ammonification in suboxic waters with denitrification
         ! -------------------------------------------------------
         zammonic = zremikc * nitrfac * doc
         denitr  = zammonic * ( 1._rk - nitrfac2 )
         denitr  = MAX(0._rk, MIN(  ( no3 - rtrn ) / rdenit, denitr ) )
         zoxyremc          = MAX(0., zammonic - denitr)
         !

         _ADD_SOURCE_(self%id_po4, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_nh4, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_no3, - denitr * rdenit)
         _ADD_SOURCE_(self%id_doc, - zolimi - denitr - zoxyremc)
         _ADD_SOURCE_(self%id_oxy, - zolimi * o2ut)
         _ADD_SOURCE_(self%id_dic, + zolimi + denitr + zoxyremc)
         _ADD_SOURCE_(self%id_tal, + rno3 * ( zolimi + zoxyremc + ( rdenit + 1._rk) * denitr ))

         _SET_DIAGNOSTIC_(self%id_remin, zolimi)
         _SET_DIAGNOSTIC_(self%id_denit, denitr * rdenit * rno3)

         zdep = MAX( hmld, heup_01) !, gdept_n)
         zdepmin = MIN( 1._rk, zdep / gdept_n )
         zdepprod = zdepmin**0.273_rk
         zdepeff  = 0.3_rk * zdepmin**0.6_rk
         tgfunc = EXP( 0.0631_rk * tem ) 


         !--------------

         ! Bacterial uptake of iron. No iron is available in DOC. So
         ! Bacteries are obliged to take up iron from the water. Some
         ! studies (especially at Papa) have shown this uptake to be significant
         ! ----------------------------------------------------------
         zbactfer = self%feratb * self%mumax0 * xstep * tgfunc * xlimbacl     & ! Jorn: Eq 63, mu now configurable (defaulting to 0.6, OA 2021-09-03), dropped multiplication with rfact2 [time step in seconds]
           &              * fer / ( self%xkferb + fer )    &
           &              * zdepeff * zdepbac                          ! Mokrane: in the original version we use biron instead of fer
         _ADD_SOURCE_(self%id_fer, - zbactfer*0.1_rk)
         _ADD_SOURCE_(self%id_sfe, + zbactfer*0.08_rk)
         _ADD_SOURCE_(self%id_bfe, + zbactfer*0.02_rk)
         _SET_DIAGNOSTIC_(self%id_febact, zbactfer * 0.1_rk)
         blim      = xlimbacl  * zdepbac / 1.e-6_rk 
         _SET_DIAGNOSTIC_(self%id_blim, blim)

      _LOOP_END_
   end subroutine

   subroutine bacteria_initialize(self, configunit)
      class (type_pisces_bacteria), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))

      call self%register_diagnostic_variable(self%id_bact, 'bact', 'mol C L-1', 'biomass proxy', source=source_do_column)

      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'euphotic depth')
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      !call self%register_dependency(self%id_hmld, standard_variables%mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_dependency(self%id_mes, 'mes', 'mol C L-1', 'mesozooplankton')
   end subroutine

   subroutine bacteria_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_bacteria), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: gdept_n, hmld, heup_01, zoo, mes
      real(rk) :: zdep, zdepbac, ztempbac, zdepmin

      _DOWNWARD_LOOP_BEGIN_
         ! Jorn from p4zrem.F90
         ! Computation of the mean phytoplankton concentration as
         ! a crude estimate of the bacterial biomass
         ! this parameterization has been deduced from a model version
         ! that was modeling explicitely bacteria
         ! -------------------------------------------------------
         _GET_(self%id_gdept_n, gdept_n)
         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_SURFACE_(self%id_heup_01, heup_01)
         _GET_(self%id_zoo, zoo)
         _GET_(self%id_mes, mes)
         zdep = MAX( hmld, heup_01 )      ! Jorn: Eq 35a
         IF( gdept_n < zdep ) THEN
            zdepbac = 0.6_rk * ( MAX(0._rk, zoo + mes ) * 1.0E6 )**0.6 * 1.E-6 !MIN( 0.7_rk * ( zoo + 2._rk* mes ), 4.e-6_rk )    ! Jorn: Eq 35b
            ztempbac   = zdepbac     ! Jorn: this saves the last value that is then used in the deep (clause below)
         ELSE
            zdepmin = MIN( 1._rk, zdep / gdept_n )
            zdepbac = zdepmin**0.683_rk * ztempbac    ! Jorn: Eq 35b
         ENDIF
         _SET_DIAGNOSTIC_(self%id_bact, zdepbac)
      _DOWNWARD_LOOP_END_
   end subroutine

end module
