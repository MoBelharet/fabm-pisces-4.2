#include "fabm_driver.h"

module pisces_zooplankton

   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_zooplankton
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_don, id_dop, id_oxy, id_fer
      type (type_state_variable_id) :: id_dia, id_dfe, id_dsi, id_dch, id_phy, id_nfe, id_nch, id_zoo, id_poc, id_sfe, id_goc, id_bfe, id_gsi
      type (type_state_variable_id) :: id_po4, id_no3, id_nh4, id_doc, id_dic, id_tal, id_poc_waste, id_pof_waste, id_pos_waste, id_cal
      type (type_state_variable_id) :: id_conspoc, id_consgoc, id_prodpoc, id_poc_waste_prod, id_prodgoc
      type (type_dependency_id)     :: id_tem, id_nitrfac, id_quotan, id_quotad, id_xfracal, id_wspoc, id_wsgoc, id_gdepw_n, id_sized, id_sizen
      !type (type_global_dependency_id)  :: id_rDttrc
      type (type_surface_dependency_id) :: id_hmld, id_heup_01
      type (type_diagnostic_variable_id) :: id_zfezoo, id_zgrazing, id_zfrac, id_pcal
      type (type_diagnostic_variable_id) :: id_zfood_diag, id_zfoodlim_diag, id_ztmp2_diag, id_ztmp1_diag, id_zgrazp_diag, id_zgrazpoc_diag,&
      & id_zcompaph_diag, id_zcompapoc_diag, id_zcompadi_diag, id_sizen_diag, id_sized_diag , id_zgraze_diag,id_zdiffdn_diag,id_ztmp3_diag, id_ztmp4_diag, &
      & id_zproport_diag, id_zgrasrat_diag, id_zgrasratn_diag, id_zepshert_diag, id_zepsherv_diag, id_quotan_diag, id_zepsherf_diag, &
      & id_zepsherq_diag,  id_zbeta_diag, id_zgrazpof_diag, id_zgrazd_diag, id_zgraztotf_diag, id_zgraztotc_diag, id_zgraztotn_diag

      real(rk) :: grazrat, logbz, resrat, xkmort, mzrat, xthresh, xkgraz, ferat, feratn, epsher, epshermin, unass, sigma, part, grazflux
      real(rk) :: xthreshdia, xthreshphy, xthreshzoo, xthreshpoc
      real(rk) :: xprefn, xprefz, xprefd, xprefc, xsizedia, xdismort, phlim, xsigma, xsigmadel
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_pisces_zooplankton), intent(inout), target :: self
      integer,                         intent(in)            :: configunit

      real(rk) :: bz

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(self%epsher, 'epsher', '1', 'maximum growth efficiency', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%epshermin, 'epshermin', '1', 'minimum growth efficiency', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%unass, 'unass', '1', 'non-assimilated fraction', default=0.3_rk, minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%sigma, 'sigma', '1', 'excretion as dissolved organic matter', default=0.6_rk, minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%grazrat, 'grazrat', 'd-1', 'maximum grazing rate', minimum=0._rk)
      call self%get_parameter(bz, 'bz', '-', 'Temperature sensitivity of grazing rate',default=1.079_rk)
      call self%get_parameter(self%logbz, 'logbz', '-', 'Temperature sensitivity of grazing rate (log rate, overwrites bz if given explicitely)', default= log(bz))
      call self%get_parameter(self%grazflux, 'grazflux', '(m mol C L-1)-1', 'flux-feeding rate', minimum=0._rk)
      call self%get_parameter(self%xkgraz, 'xkgraz', 'mol C L-1', 'half-saturation constant for grazing', default=20.e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xprefn, 'xprefn', '-', 'preference for nanophytoplankton', minimum=0._rk)
      call self%get_parameter(self%xprefd, 'xprefd', '-', 'preference for diatoms', minimum=0._rk)
      call self%get_parameter(self%xprefz, 'xprefz', '-', 'preference for microzooplankton', minimum=0._rk)
      call self%get_parameter(self%xprefc, 'xprefc', '-', 'preference for particulate organic carbon', minimum=0._rk)
      call self%get_parameter(self%xthresh, 'xthresh', 'mol C L-1', 'food threshold', default=0.3e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshdia, 'xthreshdia', 'mol C L-1', 'food threshold for diatoms', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshphy, 'xthreshphy', 'mol C L-1', 'food threshold for nanophytoplankton', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshzoo, 'xthreshzoo', 'mol C L-1', 'food threshold for microzooplankton', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%xthreshpoc, 'xthreshpoc', 'mol C L-1', 'food threshold for particulate organic carbon', default=1e-8_rk, minimum=0._rk)
      call self%get_parameter(self%mzrat, 'mzrat', '(umol C L-1)-1 d-1', 'quadratic mortality', minimum=0._rk)
      call self%get_parameter(self%resrat, 'resrat', 'd-1', 'linear mortality', minimum=0._rk)
      call self%get_parameter(self%xkmort, 'xkmort', 'mol C L-1', 'half saturation constant for mortality', default=0.1e-6_rk, minimum=0._rk) ! UPDATED
      call self%get_parameter(self%part, 'part', '1', 'fraction of calcite that does not dissolve in guts', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%ferat, 'ferat', 'mol Fe mol C-1', 'Fe / C ratio', default=10.e-6_rk, minimum=0._rk)
      call self%get_parameter(self%feratn, 'feratn', 'mol Fe mol C-1', 'Fe / C ratio of microzooplankton', default=10.e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xsizedia, 'xsizedia', 'mol C L-1', 'maximum accessible diatom biomass (above this threshold cells are too large)', default=1e-6_rk, minimum=0._rk)
      call self%get_parameter(self%xdismort, 'xdismort', '1', 'fraction of quadratic mortality directed to nutrient pools', default=( 1._rk - self%epsher - self%unass ) /( 1._rk - self%epsher ), minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%phlim, 'phlim', '1', 'relative grazing on nanophytoplankton if cells are small', minimum=0._rk, maximum=1._rk)
      call self%get_parameter(self%xsigma, 'xsigma', '1', 'Width of the grazing window' , default=0.5_rk , minimum=0._rk) 
      call self%get_parameter(self%xsigmadel, 'xsigmadel', '1', 'Maximum additional width of the grazing window', default=1._rk , minimum=0._rk)

      call self%register_state_variable(self%id_c, 'c', 'mol C L-1', 'carbon', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rno3 * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=po4r * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron, self%id_c, scale_factor=self%ferat * 1e9_rk)

      call self%register_diagnostic_variable(self%id_zfezoo, 'fezoo', 'nmol Fe m-3 s-1', 'iron recycling rate')
      call self%register_diagnostic_variable(self%id_zgrazing, 'graz', 'mol C m-3 s-1', 'grazing')
      call self%register_diagnostic_variable(self%id_zfrac, 'zfrac', 'mol C m-3 s-1', 'fractionation of large POM')
      call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mol m-3 s-1', 'calcite production')
      call self%add_to_aggregate_variable(calcite_production, self%id_pcal)

      call self%register_diagnostic_variable(self%id_zfood_diag, 'zfood_diag', '-', 'zfood diagnostic')
      call self%register_diagnostic_variable(self%id_zfoodlim_diag, 'zfoodlim_diag', '-', 'zfoodlim diagnostic')
      call self%register_diagnostic_variable(self%id_ztmp2_diag, 'ztmp2_diag', '-', 'ztmp2 diagnostic')
      call self%register_diagnostic_variable(self%id_ztmp1_diag, 'ztmp1_diag', '-', 'ztmp1 diagnostic')
      call self%register_diagnostic_variable(self%id_ztmp3_diag, 'ztmp3_diag', '-', 'ztmp3 diagnostic')
      call self%register_diagnostic_variable(self%id_ztmp4_diag, 'ztmp4_diag', '-', 'ztmp4 diagnostic')
      call self%register_diagnostic_variable(self%id_zgrazp_diag, 'zgrazp_diag', '-', 'zgrazp diagnostic')
      call self%register_diagnostic_variable(self%id_zgrazpoc_diag, 'zgrazpoc_diag', '-', 'zgrazpoc diagnostic')
      call self%register_diagnostic_variable(self%id_zcompaph_diag, 'zcompaph_diag', '-', 'zcompaph diagnostic')
      call self%register_diagnostic_variable(self%id_zcompapoc_diag, 'zcompapoc_diag', '-', 'zcompapoc diagnostic')
      call self%register_diagnostic_variable(self%id_zcompadi_diag, 'zcompadi_diag', '-', 'zcompadi diagnostic')
      call self%register_diagnostic_variable(self%id_sizen_diag, 'sizen_diag' , '-', 'diagnostic of sizen' )
      call self%register_diagnostic_variable(self%id_sized_diag, 'sized_diag' , '-', 'diagnostic of sized' )
      call self%register_diagnostic_variable(self%id_zgraze_diag, 'zgraze_diag' , '-', 'diagnostic of zgraze' )
      call self%register_diagnostic_variable(self%id_zdiffdn_diag, 'zdiffdn_diag', '-', 'diagnostic of zdiffdn' )
      call self%register_diagnostic_variable(self%id_zproport_diag, 'zproport_diag', '-','diagnostic of zproport')
      call self%register_diagnostic_variable(self%id_zgrasrat_diag, 'zgrasrat_diag', '-','diagnostic of zgrasrat')
      call self%register_diagnostic_variable(self%id_zgrasratn_diag, 'zgrasratn_diag', '-','diagnostic of zgrasratn')
      call self%register_diagnostic_variable(self%id_zepshert_diag, 'zepshert_diag', '-','diagnostic of zepshert')
      call self%register_diagnostic_variable(self%id_zepsherv_diag, 'zepsherv_diag', '-','diagnostic of zepsherv')
      call self%register_diagnostic_variable(self%id_quotan_diag, 'quotan_diag' , '-','diagnostic of quotan')
      call self%register_diagnostic_variable(self%id_zepsherf_diag, 'zepsherf_diag', '-','diagnostic of zepsherf')
      call self%register_diagnostic_variable(self%id_zepsherq_diag, 'zepsherq_diag', '-','diagnostic of zepsherq')
      call self%register_diagnostic_variable(self%id_zbeta_diag, 'zbeta_diag', '-','diagnostic of zbeta')
      call self%register_diagnostic_variable(self%id_zgrazpof_diag, 'zgrazpof_diag', '-','diagnostic of zgrazpof')
      call self%register_diagnostic_variable(self%id_zgrazd_diag, 'zgrazd_diag', '-','diagnostic of zgrazd')
      call self%register_diagnostic_variable(self%id_zgraztotf_diag, 'zgraztotf_diag', '-', 'diagnostic of zgraztotf' )
      call self%register_diagnostic_variable(self%id_zgraztotc_diag, 'zgraztotc_diag', '-', 'diagnostic of zgraztotc' )
      call self%register_diagnostic_variable(self%id_zgraztotn_diag, 'zgraztotn_diag', '-', 'diagnostic of zgraztotn' )

      call self%register_state_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_state_dependency(self%id_nh4, 'nh4', 'mol C L-1', 'ammonium')
      call self%register_state_dependency(self%id_po4, 'po4', 'mol C L-1', 'phosphate')
      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_dic, standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)
      call self%register_state_dependency(self%id_tal, standard_variables%alkalinity_expressed_as_mole_equivalent)
      call self%register_state_dependency(self%id_poc_waste, 'poc_waste', 'mol C L-1', 'particulate carbon waste')
      call self%register_state_dependency(self%id_pof_waste, 'pof_waste', 'mol Fe L-1', 'particulate iron waste')
      call self%register_state_dependency(self%id_pos_waste, 'pos_waste', 'mol Si L-1', 'particulate silicon waste')
      call self%register_state_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')
      call self%register_state_dependency(self%id_poc_waste_prod, 'poc_waste_prod', 'mol C L-1', 'produced particulate carbon waste')
      call self%request_coupling_to_model(self%id_poc_waste, 'waste', 'c')
      call self%request_coupling_to_model(self%id_pof_waste, 'waste', 'fe')
      call self%request_coupling_to_model(self%id_pos_waste, 'waste', 'si')
      call self%request_coupling_to_model(self%id_poc_waste_prod, 'waste', 'prod')

      ! Prey
      call self%register_state_dependency(self%id_dia, 'diac', 'mol C L-1', 'diatoms')
      call self%register_state_dependency(self%id_dfe, 'dfe', 'mol Fe L-1', 'diatom iron')
      call self%register_state_dependency(self%id_dsi, 'dsi', 'mol Si L-1', 'diatom silicon')
      call self%register_state_dependency(self%id_dch, 'dch', 'g L-1', 'diatom chlorophyll')
      call self%register_state_dependency(self%id_phy, 'phyc', 'mol C L-1', 'nanophytoplankton')
      call self%register_state_dependency(self%id_nfe, 'nfe', 'mol Fe L-1', 'nanophytoplankton iron')
      call self%register_state_dependency(self%id_nch, 'nch', 'g L-1', 'nanophytoplankton chlorophyll')
      call self%register_state_dependency(self%id_zoo, 'zoo', 'mol C L-1', 'microzooplankton')
      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_state_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'large particulate organic silicon')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')
      call self%register_state_dependency(self%id_prodpoc, 'prodpoc', 'mol C L-1', 'produced small particulate organic carbon')
      call self%register_state_dependency(self%id_prodgoc, 'prodgoc', 'mol C L-1', 'produced big particulate organic carbon') !Mokrane
      call self%register_state_dependency(self%id_conspoc, 'conspoc', 'mol C L-1', 'consumed small particulate organic carbon')
      call self%register_state_dependency(self%id_consgoc, 'consgoc', 'mol C L-1', 'consumed large particulate organic carbon')
      call self%register_dependency(self%id_quotan, 'quotan', '-', 'proxy for N quota of nanophytoplankton')
      call self%register_dependency(self%id_quotad, 'quotad', '-', 'proxy for N quota of diatoms')
      call self%register_dependency(self%id_xfracal, 'xfracal', '1', 'calcifying fraction of nanophytoplankton')
      call self%register_dependency(self%id_wspoc, 'wspoc', 'm d-1', 'sinking velocity of small particulate organic matter')
      call self%register_dependency(self%id_wsgoc, 'wsgoc', 'm d-1', 'sinking velocity of large particulate organic matter')
      call self%register_dependency(self%id_sizen,'sizen','-','Mean relative size of nanophytoplankton')
      call self%register_dependency(self%id_sized,'sized','-','Mean relative size of diatoms')
      call self%register_dependency(self%id_hmld, mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      !call self%register_dependency(self%id_hmld, standard_variables%mixed_layer_thickness_defined_by_vertical_tracer_diffusivity)
      call self%register_dependency(self%id_heup_01, 'heup_01', 'm', 'euphotic layer depth (PAR > 0.5 W m-2)')
      call self%register_dependency(self%id_gdepw_n, standard_variables%depth)
      !call self%register_dependency(self%id_rDttrc, standard_variables%physical_time_step)

      ! Allow wholesale coupling of prey, which then automatically sets up the constituent coupling
      call self%request_coupling_to_model(self%id_dia, 'dia', 'c')
      call self%request_coupling_to_model(self%id_dfe, 'dia', 'fe')
      call self%request_coupling_to_model(self%id_dch, 'dia', 'ch')
      call self%request_coupling_to_model(self%id_dsi, 'dia', 'si')
      call self%request_coupling_to_model(self%id_quotad, 'dia', 'quota')
      call self%request_coupling_to_model(self%id_sized, 'dia', 'sizep')
      call self%request_coupling_to_model(self%id_phy, 'phy', 'c')
      call self%request_coupling_to_model(self%id_nfe, 'phy', 'fe')
      call self%request_coupling_to_model(self%id_nch, 'phy', 'ch')
      call self%request_coupling_to_model(self%id_quotan, 'phy', 'quota')
      call self%request_coupling_to_model(self%id_sizen, 'phy', 'sizep')
      call self%request_coupling_to_model(self%id_xfracal, 'phy', 'xfracal')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_wspoc, 'pom', 'ws')
      call self%request_coupling_to_model(self%id_prodpoc, 'pom', 'prod')
      call self%request_coupling_to_model(self%id_conspoc, 'pom', 'cons')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_gsi, 'gom', 'si')
      call self%request_coupling_to_model(self%id_cal, 'gom', 'cal')
      call self%request_coupling_to_model(self%id_wsgoc, 'gom', 'ws')
      call self%request_coupling_to_model(self%id_prodgoc, 'gom', 'prod') ! Mokrane
      call self%request_coupling_to_model(self%id_consgoc, 'gom', 'cons')

      call self%register_dependency(self%id_tem, standard_variables%temperature)
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: c, tem, nitrfac, rfact2
      real(rk) :: tgfunc2, zcompa, zfact, zrespz, ztortz !, rDt_trc
      real(rk) :: zcompadi, zcompaph, zcompapoc, zcompaz
      real(rk) :: zfood, zfoodlim, zdenom, zdenom2, zgraze
      real(rk) :: dia, dfe, dsi, dch, phy, nfe, nch, zoo, poc, sfe, goc, bfe, gsi, quotan, quotad, xfracal, wspoc, wsgoc
      real(rk) :: zgrazp, zgrazpoc, zgrazd, zgrazz, zgrazpf, zgrazpof, zgrazsf, zgraztotc, zgraztotf, zgraztotn
      real(rk) :: zgrazffeg, zgrazfffg, zgrazffep, zgrazfffp, zproport, zratio, zratio2, zfrac, zfracfe
      real(rk) :: zgrasrat, zgrasratn, zepshert, zbeta, zepsherf, zepsherq, zepsherv, zgrafer, zgrarem, zgrarsig, zmortz, zmortzgoc
      real(rk) :: zprcaca, zfracal, zgrazcal, cal, zsigma, zdiffdn
      real(rk) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmptot, sizen, sized, gdepw_n, hmld, heup_01, zgrapoc, zgrapof

      _LOOP_BEGIN_
         _GET_(self%id_c, c)
         _GET_(self%id_tem, tem)
         _GET_(self%id_nitrfac, nitrfac)
         _GET_(self%id_sizen, sizen)
         _GET_(self%id_sized, sized)

         _SET_DIAGNOSTIC_(self%id_sizen_diag, sizen)
         _SET_DIAGNOSTIC_(self%id_sized_diag, sized)

         !_GET_GLOBAL_(self%id_rDttrc, rDt_trc)

         rfact2 = rDt_trc /1._rk       !

          tgfunc2 = EXP( 0.0761_rk  * tem )        ! Jorn: from p4zint.F90, equivalent to Eq 25b as NB EXP(0.07608) = 1.079 = b_Z
!         tgfunc2 = EXP( self%logbz  * tem )         ! AC: set log(bz) as parameter.

         zcompa = MAX( ( c - 1.e-9_rk ), 0.e0_rk )
         zfact   = xstep * tgfunc2 * zcompa

         zproport  = MIN(1._rk, exp(-1.1_rk * MAX(0._rk, self%phlim * ( sized - 1.8_rk ))**0.8 )) ! for meso : zproport = 1 , 

         !   Michaelis-Menten mortality ates of microzooplankton - Jorn: last (4th) term in Eq 24
         !   -----------------------------------------------------
         zrespz = self%resrat * zfact * ( c / ( self%xkmort + c )  &
         &        + 3._rk * nitrfac )

         !   Zooplankton mortality. A square function has been selected with
         !   no real reason except that it seems to be more stable and may mimic predation.
         !   Jorn: 3rd term in Eq 24, except that a (1._rk - nitrfac) factor has been added here, and zcompaz is used to introduce a threshold
         !   ------------------------------------------------------------------------------
         ztortz = self%mzrat * 1.e6_rk * zfact * c * (1._rk - nitrfac)

         ! Prey
         _GET_(self%id_dia, dia)
         _GET_(self%id_dfe, dfe)
         _GET_(self%id_dsi, dsi)
         _GET_(self%id_dch, dch)
         _GET_(self%id_phy, phy)
         _GET_(self%id_nfe, nfe)
         _GET_(self%id_nch, nch)
         _GET_(self%id_zoo, zoo)
         _GET_(self%id_poc, poc)
         _GET_(self%id_sfe, sfe)
         _GET_(self%id_goc, goc)
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_bfe, bfe)
         _GET_(self%id_cal, cal)
         _GET_(self%id_wspoc, wspoc)
         _GET_(self%id_wsgoc, wsgoc)
         _GET_(self%id_quotan, quotan)
         _GET_(self%id_quotad, quotad)
         _GET_(self%id_xfracal, xfracal)
       
         zcompadi  = zproport * MAX( ( dia - self%xthreshdia ), 0._rk )  ! Jorn: xsizedia = 0 for mesozoo
         zcompaph  = MAX( ( phy - self%xthreshphy ), 0._rk )
         zcompapoc = MAX( ( poc - self%xthreshpoc ), 0._rk )
         zcompaz   = MAX( ( zoo - self%xthreshzoo ), 0._rk )

         _SET_DIAGNOSTIC_(self%id_zcompaph_diag, zcompaph)
         _SET_DIAGNOSTIC_(self%id_zcompapoc_diag, zcompapoc)
         _SET_DIAGNOSTIC_(self%id_zcompadi_diag, zcompadi)

         ! Jorn - mesozo only:
         ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
         ! it is to predation by mesozooplankton
         ! -------------------------------------------------------------------------------

         ! Grazing
         ! ----------------------------------
         zfood     = self%xprefn * zcompaph + self%xprefc * zcompapoc + self%xprefd * zcompadi + self%xprefz * zcompaz   ! Jorn: 1st term in Eq 26a
         _SET_DIAGNOSTIC_(self%id_zfood_diag, zfood)
         zfoodlim  = MAX( 0._rk , zfood - MIN(self%xthresh, 0.5_rk * zfood) )                                               ! Jorn: 2nd term in Eq 26a
         _SET_DIAGNOSTIC_(self%id_zfoodlim_diag, zfoodlim)
         zdenom    = zfoodlim / ( self%xkgraz + zfoodlim )
         zdenom2   = zdenom / ( zfood + rtrn )
         zgraze    = self%grazrat * xstep * tgfunc2 * c * (1. - nitrfac)!     ! Jorn: compared to paper (Eq 26a), (1._rk - nitrfac) factor seems to have been added
          

         _SET_DIAGNOSTIC_(self%id_zgraze_diag, zgraze/xstep)

         zsigma = 1.0_rk - zdenom**2/(0.05_rk**2 + zdenom**2)
         zsigma = self%xsigma + self%xsigmadel * zsigma

         zdiffdn = exp( -ABS(log(1.67_rk * sizen / (5._rk * sized + rtrn )) )**2 / zsigma**2)

         _SET_DIAGNOSTIC_(self%id_zdiffdn_diag, zdiffdn)

         ztmp1 = self%xprefn * zcompaph * ( zcompaph + zdiffdn * zcompadi ) / ( 1._rk + zdiffdn ) !
         ztmp2 = self%xprefd * zcompadi * ( zcompadi + zdiffdn * zcompaph ) / ( 1._rk + zdiffdn ) !
         ztmp3 = self%xprefc * zcompapoc * zcompapoc !
         ztmp4 = self%xprefz * zcompaz * zcompaz !

         _SET_DIAGNOSTIC_(self%id_ztmp1_diag, ztmp1)
         _SET_DIAGNOSTIC_(self%id_ztmp2_diag, ztmp2)
         _SET_DIAGNOSTIC_(self%id_ztmp3_diag, ztmp3)
         _SET_DIAGNOSTIC_(self%id_ztmp4_diag, ztmp4)

         ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + rtrn ! 
         ztmp1 = ztmp1 / ztmptot
         ztmp2 = ztmp2 / ztmptot
         ztmp3 = ztmp3 / ztmptot
         ztmp4 = ztmp4 / ztmptot

         zgrazd    = zgraze  * ztmp2 * zdenom  ! diatoms   ! A déclarer
         zgrazp    = zgraze  * ztmp1 * zdenom  ! nanophytoplankton ! A déclarer
         zgrazpoc  = zgraze  * ztmp3 * zdenom  ! small POC
         zgrazz    = zgraze  * ztmp4 * zdenom  ! microzooplankton


         zgrazpf   = zgrazp   * nfe / ( phy + rtrn)  !a
         zgrazsf   = zgrazd   * dfe / ( dia + rtrn)  ! 
         zgrazpof  = zgrazpoc * sfe / ( poc + rtrn)

         _SET_DIAGNOSTIC_(self%id_zgrazpof_diag, zgrazpof/xstep)

         ! Flux feeding on GOC
         ! ----------------------------------
         zgrazffeg = self%grazflux  * xstep * wsgoc      &  ! Jorn: Eq 29b but with added factor (1._rk - nitrfac). NB sinking velocity wsgoc should be m d-1!
         &           * tgfunc2 * goc * c &                  ! for microzoo : zgrazffeg = 0
         &           * (1._rk - nitrfac)
         zgrazfffg = zgrazffeg * bfe / (goc + rtrn)
         zgrazffep = self%grazflux  * xstep *  wspoc     &  ! Jorn: Eq 29a but with added factor (1._rk - nitrfac). NB sinking velocity wsgoc should be m d-1!
         &           * tgfunc2 * poc * c &                  ! for microzoo : zgrazffep = 0
         &           * (1._rk - nitrfac)
         zgrazfffp = zgrazffep * sfe / (poc + rtrn)
         !
         zgraztotc = zgrazp + zgrazd  + zgrazpoc & ! Jorn: this seems to be potential ingestion, since zgrazffep and zgrazffeg will later be scaled with proportion of filter feeders, zproport
                 &   + zgrazz + zgrazffep + zgrazffeg   ! Mokrane: this term only for meso (zgrazz = zgrazffep = zgrazffeg = 0 for micro)
      
         !********** Mokrane: This part concerns only Mesozoo, nothing will be changed for microzoo ************************************************************

         _GET_SURFACE_(self%id_hmld, hmld)
         _GET_SURFACE_(self%id_heup_01, heup_01)
         _GET_(self%id_gdepw_n, gdepw_n)
         zproport  = 0._rk
         IF( gdepw_n > MAX(hmld , heup_01 ) ) THEN
              zproport  = (zgrazffep + zgrazffeg)/(rtrn/rfact2 + zgraztotc) ! for microzoo: zproport = 0 because zgrazffep=0 and zgrazffeg=0
         ENDIF
         _SET_DIAGNOSTIC_(self%id_zproport_diag, zproport)


         zratio    = gsi / ( goc + rtrn )
         zratio2   = zratio * zratio
         zfrac     = zproport * self%grazflux  * xstep * wsgoc      &
         &          * goc * c          &
         &          * ( 0.4_rk + 3.6_rk * zratio2 / ( 1._rk**2 + zratio2 ) )  ! (zfrac=0 for microzoo)
         zfracfe   = zfrac * bfe / (goc + rtrn)                               ! (zfracfe=0 for mesozoo)


         zgrazffep = zproport * zgrazffep ! for micro: zgrazffep=0
         zgrazffeg = zproport * zgrazffeg ! for micro: zgrazffeg =0
         zgrazfffp = zproport * zgrazfffp ! for micro: zgrazfffp =0
         zgrazfffg = zproport * zgrazfffg ! for micro: zgrazfffg =0

         ! Mokrane: zproport = zproport for mesozoo and zproport = 0 for microzoo

         zgrazd    = (1.0 - zproport) * zgrazd ! for micro: zgrazd = zgrazd because zproport = 0
         zgrazp    = (1.0 - zproport) * zgrazp
         zgrazz    = (1.0 - zproport) * zgrazz
         zgrazpoc  = (1.0 - zproport) * zgrazpoc
         zgrazsf   = (1.0 - zproport) * zgrazsf
         zgrazpf   = (1.0 - zproport) * zgrazpf
         zgrazpof  = (1.0 - zproport) * zgrazpof


         !*********************************************************************************************************************************************************

         zgraztotc = zgrazp  + zgrazd  + zgrazpoc &
                 &   + zgrazz + zgrazffep + zgrazffeg    ! Mokrane: this term = 0 for micro
         zgraztotf = zgrazpf + zgrazsf + zgrazpof &
                 &   + zgrazz * self%feratn   + zgrazfffp + zgrazfffg    ! Mokrane: this term = 0 for micro
         zgraztotn = zgrazp * quotan + zgrazd * quotad + zgrazpoc &
                 &   + zgrazz  + zgrazffep + zgrazffeg ! Mokrane: this term = 0 for micro

         _SET_DIAGNOSTIC_(self%id_quotan_diag, quotan)

         ! Grazing by microzooplankton
         _SET_DIAGNOSTIC_(self%id_zgrazing, zgraztotc * 1e3_rk)

         ! Microzooplankton efficiency.
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and
         ! Fulton, 2012)
         ! -----------------------------------------------------------------------------
         zgrasrat  = ( zgraztotf + rtrn/rfact2 ) / ( zgraztotc + rtrn/rfact2  )  ! Jorn: Fe : C ratio in ingested prey
         zgrasratn = ( zgraztotn + rtrn/rfact2 ) / ( zgraztotc + rtrn/rfact2  )  ! Jorn: N : C ratio in ingested prey, but N already expressed in C units
         zepshert  =  MIN( 1._rk, zgrasratn, zgrasrat / self%ferat)  ! Jorn: Eq 27a, maximum rate of biomass production, derived from incoming C, N, Fe
         zbeta     = MAX(0._rk, (self%epsher - self%epshermin) )
         zepsherf  = self%epshermin + zbeta / ( 1.0_rk + 0.04E6_rk * 12._rk * zfood * zbeta )
         zepsherq  = 0.5_rk + (1.0_rk - 0.5_rk) * zepshert * ( 1.0_rk + 1.0_rk ) / ( zepshert + 1.0_rk )
         zepsherv  = zepsherf * zepshert * zepsherq   ! Jorn: gross growth efficiency

         _SET_DIAGNOSTIC_(self%id_zgraztotf_diag , zgraztotf/xstep)
         _SET_DIAGNOSTIC_(self%id_zgraztotc_diag , zgraztotc/xstep)
         _SET_DIAGNOSTIC_(self%id_zgraztotn_diag , zgraztotn/xstep)
         _SET_DIAGNOSTIC_(self%id_zgrasrat_diag, zgrasrat)
         _SET_DIAGNOSTIC_(self%id_zgrasratn_diag, zgrasratn)
         _SET_DIAGNOSTIC_(self%id_zepshert_diag, zepshert)
         _SET_DIAGNOSTIC_(self%id_zepsherv_diag, zepsherv)
         _SET_DIAGNOSTIC_(self%id_zepsherf_diag, zepsherf)
         _SET_DIAGNOSTIC_(self%id_zepsherq_diag, zepsherq)
         _SET_DIAGNOSTIC_(self%id_zbeta_diag, zbeta)

!****************************************************************************************************************************************
         zmortz = ztortz + zrespz
         zmortzgoc = (1._rk - self%xdismort) * ztortz + zrespz

         _ADD_SOURCE_(self%id_c, - zmortz + zepsherv * zgraztotc )
         _ADD_SOURCE_(self%id_zoo, - zgrazz)
         !---------------------------------------------------------
         _ADD_SOURCE_(self%id_phy, - zgrazp)
         _ADD_SOURCE_(self%id_dia, - zgrazd)
         _SET_DIAGNOSTIC_(self%id_zgrazd_diag, zgrazd/xstep)
         _SET_DIAGNOSTIC_(self%id_zgrazp_diag, zgrazp/xstep)
         _SET_DIAGNOSTIC_(self%id_zgrazpoc_diag, zgrazpoc/xstep)
         _ADD_SOURCE_(self%id_nch, - zgrazp  * nch / (phy + rtrn))
         _ADD_SOURCE_(self%id_dch, - zgrazd * dch / (dia + rtrn))
         _ADD_SOURCE_(self%id_dsi, - zgrazd * dsi / (dia + rtrn))
         !----------Mokrane --------

         _ADD_SOURCE_(self%id_gsi, + zgrazd * dsi / (dia + rtrn))

         !-------------------------------
         !zgrabsi       = zgrazd * dsi / ( dia + rtrn ) ! zgrabsi to be declared
         _ADD_SOURCE_(self%id_nfe, - zgrazpf)
         _ADD_SOURCE_(self%id_dfe, - zgrazsf)
         _ADD_SOURCE_(self%id_poc, - zgrazpoc - zgrazffep + zfrac)
         !------------------------------------------------------------------
         _ADD_SOURCE_(self%id_conspoc, - zgrazpoc - zgrazffep)
         _ADD_SOURCE_(self%id_prodpoc, + zfrac)
         _ADD_SOURCE_(self%id_goc,            - zgrazffeg - zfrac )
         _ADD_SOURCE_(self%id_consgoc,        - zgrazffeg - zfrac)
         _ADD_SOURCE_(self%id_sfe, - zgrazpof - zgrazfffp + zfracfe)
         _ADD_SOURCE_(self%id_bfe,            - zgrazfffg - zfracfe)
         zfracal = cal / (goc + rtrn )
         zgrazcal = zgrazffeg * (1._rk - self%part) * zfracal ! for microzoo: zgrazcal = 0 because zgrazffeg = 0
         zprcaca = xfracal * zgrazp
         zprcaca = self%part * zprcaca

         zgrafer   = zgraztotc * MAX( 0._rk , ( 1._rk - self%unass ) * zgrasrat   - self%ferat * zepsherv ) &
!         zgrafer   = ( 1._rk - self%unass ) * zgraztotf - self%ferat * zepsherv * zgraztotc &  ! Jorn: total dissolved Fe waste! (organic + inorganic). TODO: revert to original eq above? non-conservative!
         &         + self%ferat * self%xdismort * ztortz * (1._rk - self%phlim)! Jorn: Eq30b. this line for mesozoo only , ztortz is quadratic mortality

         zgrarem   = zgraztotc * ( 1._rk - zepsherv - self%unass ) &   ! Jorn: total dissolved C/N/P waste (organic + inorganic)
                 &         + self%xdismort * ztortz * (1._rk - self%phlim) ! Jorn: Eq30b. this line for mesozoo only (1 - phlim = 0),

         zgrapoc  = zgraztotc * self%unass &
                 &  + zmortzgoc * (1._rk - self%phlim)  ! Mokrane: This line for mesozoo only

         zgrapof = (zgraztotf * self%unass + self%ferat * zmortzgoc)  ! Mokrane: This line for mesozoo only

         zgrarsig  = zgrarem * self%sigma


         _ADD_SOURCE_(self%id_po4, zgrarsig)
         _ADD_SOURCE_(self%id_nh4, zgrarsig)
         _ADD_SOURCE_(self%id_doc, zgrarem - zgrarsig)

         _ADD_SOURCE_(self%id_oxy, - o2ut * zgrarsig)
         _ADD_SOURCE_(self%id_fer, + zgrafer)

         _SET_DIAGNOSTIC_(self%id_zfezoo, zgrafer * 1e12_rk)

         _ADD_SOURCE_(self%id_poc, + (zgrapoc + zmortz) * self%phlim)
         _ADD_SOURCE_(self%id_goc, + zgrapoc * (1._rk - self%phlim))

         _ADD_SOURCE_(self%id_sfe, + (zgraztotf * self%unass + self%ferat * zmortz) * self%phlim)
         _ADD_SOURCE_(self%id_bfe, + zgrapof *  (1._rk-self%phlim) ) ! only for mesozoo
           
         _ADD_SOURCE_(self%id_dic, + zgrazcal - zprcaca + zgrarsig)
         _ADD_SOURCE_(self%id_tal, rno3 * zgrarsig  + 2._rk * (zgrazcal - zprcaca))
         _ADD_SOURCE_(self%id_cal, - zgrazcal + zprcaca)

!*****************************************************************************************************************


         ! Particulate waste from feeding and mortality
         !_ADD_SOURCE_(self%id_poc_waste, + zgrapoc)

         !_ADD_SOURCE_(self%id_pos_waste, + zgrazd * dsi / (dia + rtrn))
         !_ADD_SOURCE_(self%id_pof_waste, + zgrapof ) ! * (1._rk-self%phlim) ) ! Mokrane

         !_ADD_SOURCE_(self%id_poc_waste_prod, + zgraztotc * self%unass + zmortzgoc)
         !_ADD_SOURCE_(self%id_poc_waste_prod,+ zgrapoc)
         _ADD_SOURCE_(self%id_prodpoc , + (zgrapoc + zmortz) * self%phlim)
         _ADD_SOURCE_(self%id_prodgoc , +  zgrapoc * (1._rk - self%phlim))

         _SET_DIAGNOSTIC_(self%id_zfrac, + zfrac * 1e3_rk)
         !
         _SET_DIAGNOSTIC_(self%id_pcal, zprcaca * 1e3_rk)
      _LOOP_END_
   end subroutine

end module
