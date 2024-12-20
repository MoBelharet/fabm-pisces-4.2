#include "fabm_driver.h"

module pisces_iron
   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_iron
      type (type_state_variable_id) :: id_fer, id_sfe, id_bfe
      type (type_dependency_id) :: id_tempis, id_salinprac, id_xdiss, id_doc, id_poc, id_goc, id_cal, id_gsi, id_hi, id_oxy, id_etot, id_gdept_n
      type (type_dependency_id) :: id_zdust, id_etot_ndcy,id_chemo2,id_nitrfac, id_no3, id_consfe3_sum
      type (type_surface_dependency_id) :: id_gphit, id_fr_i
      type (type_diagnostic_variable_id) :: id_scav, id_coll, id_Fe3, id_FeL1, id_zTL1,id_xfecolagg, id_plig, id_zfeprecip, id_xcoagfe
      type (type_diagnostic_variable_id) :: id_zkeq_diag, id_zklight_diag, id_zconsfe_diag, id_za1_diag, id_ztfe_diag, id_zdust_diag
      real(rk) :: ligand, xlam1, xlamdust, kfep, wdust, light, scaveff
   contains
      procedure :: initialize
      procedure :: do
   end type

   logical, parameter :: ln_ligvar = .false.
   logical, parameter :: ln_ligand = .false.

contains

   subroutine initialize(self, configunit)
      class (type_pisces_iron), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(self%ligand, 'ligand', 'mol L-1', 'total concentration of iron ligands', default=1.E-9_rk)
      call self%get_parameter(self%xlam1, 'xlam1', 'd-1 umol-1 L', 'scavenging rate', default=0.02_rk)
      call self%get_parameter(self%xlamdust, 'xlamdust', 'd-1 mg-1 L', 'scavenging rate of dust', default=150.0_rk)
      call self%get_parameter(self%kfep, 'kfep', 'd-1', 'nanoparticle formation rate constant', default=0.01_rk)
      call self%get_parameter(self%light, 'light', 'W m-2', 'light limitation parameter for photolysis', default=30._rk)
      call self%get_parameter(self%scaveff, 'scaveff' , '1' , 'Fraction of scavenged Fe that goes to POFe', default=1._rk)

      call self%register_state_dependency(self%id_fer, 'fer', 'mol Fe L-1', 'iron')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol Fe L-1', 'small particulate organic iron')
      call self%register_state_dependency(self%id_bfe, 'bfe', 'mol Fe L-1', 'large particulate organic iron')

      call self%register_diagnostic_variable(self%id_scav, 'scav', 'mol Fe L-1 s-1', 'scavenging')
      call self%register_diagnostic_variable(self%id_coll, 'coll', 'mol Fe L-1 s-1', 'colloidal pumping of FeL')
      call self%register_diagnostic_variable(self%id_Fe3, 'Fe3', 'nmol Fe L-1', 'iron III concentration')
      call self%register_diagnostic_variable(self%id_FeL1, 'FeL1', 'nmol L-1', 'complexed iron concentration with L1')
      call self%register_diagnostic_variable(self%id_zTL1, 'zTL1', 'nmol L-1', 'total L1 concentration')
      call self%register_diagnostic_variable(self%id_xfecolagg,'xfecolagg','1' ,'Refractory diagnostic concentration of ligands')
      call self%register_diagnostic_variable(self%id_plig, 'plig', '1','Fraction of ligands')
      call self%register_diagnostic_variable(self%id_zfeprecip, 'zfeprecip','mmolFe L-1 s-1' ,'Precipitation of Fe')
      call self%register_diagnostic_variable(self%id_xcoagfe,'xcoagfe','mol C L-1','coagulation of colloidal iron')
      call self%register_diagnostic_variable(self%id_zkeq_diag, 'zkeq', '-', 'zkeq in iron.F90') 
      call self%register_diagnostic_variable(self%id_zklight_diag, 'zklight','-', 'zklight in iron.F90')
      call self%register_diagnostic_variable(self%id_zconsfe_diag, 'zconsfe', '-', 'zconsfe in iron.F90')
      call self%register_diagnostic_variable(self%id_za1_diag, 'za1', '-', 'za1 in iron.F90')
      call self%register_diagnostic_variable(self%id_ztfe_diag, 'ztfe', '-', 'ztfe in iron.F90') 
      call self%register_diagnostic_variable(self%id_zdust_diag, 'zdust_diag','mol m-3 s-1', 'diagnostic of zdust')

      call self%register_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_dependency(self%id_goc, 'goc', 'mol C L-1', 'large particulate organic carbon')
      call self%register_dependency(self%id_cal, 'cal', 'mol C L-1', 'calcite')
      call self%register_dependency(self%id_gsi, 'gsi', 'mol Si L-1', 'large particulate organic silicon')
      call self%register_dependency(self%id_no3, 'no3', 'mol C L-1', 'nitrate')
      call self%register_dependency(self%id_nitrfac, 'nitrfac', '1', 'denitrification factor computed from O2 levels')

      call self%request_coupling_to_model(self%id_doc, 'dom', 'c')
      call self%request_coupling_to_model(self%id_poc, 'pom', 'c')
      call self%request_coupling_to_model(self%id_goc, 'gom', 'c')
      call self%request_coupling_to_model(self%id_sfe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_bfe, 'gom', 'fe')
      call self%request_coupling_to_model(self%id_cal, 'gom', 'cal')
      call self%request_coupling_to_model(self%id_gsi, 'gom', 'si')
      call self%register_dependency(self%id_xdiss, shear_rate)
      call self%register_dependency(self%id_hi, 'hi', 'mol L-1', 'hydrogen ion concentration')
      call self%register_dependency(self%id_oxy, 'oxy', 'mol O2 L-1', 'oxygen')
      call self%register_dependency(self%id_etot, 'etot', 'W m-2', 'instantaneous PAR')
      call self%register_dependency(self%id_gphit, standard_variables%latitude)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_tempis, standard_variables%temperature) ! TODO should be in-situ temperature (as opposed to conservative/potential)
      call self%register_dependency(self%id_salinprac, standard_variables%practical_salinity)
      call self%register_dependency(self%id_zdust, 'pdust', 'g m-3', 'dust concentration')
      call self%register_dependency(self%id_etot_ndcy, 'etot_ndcy', 'W m-2', 'daily mean PAR')
      call self%register_dependency(self%id_fr_i, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_chemo2, 'chemo2', 'mol O2 (L atm)-1', 'solubility')
      call self%register_dependency(self%id_consfe3_sum, consfe3_sum)

   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_pisces_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fer, doc, poc, goc, cal, gsi, hi, oxy, etot, xdiss, gphit, tempis, salinprac, gdept_n, no3
      real(rk) :: ztkel, zsal, zis, fekeq, ztkel1, fesol(5)
      real(rk) :: ztotlig, zTL1, zkeq, zfesatur, ztfe, zFe3, zFeL1, zdust, zhplus, fe3sol, zfeequi, zfecoll, precip, ztrc, precipno3, zfeprecip
      real(rk) :: zxlam, zlam1a, zlam1b, zscave, zdenom1, zdenom2,  zcoag, zaggdfea, zaggdfeb
      real(rk) :: etot_ndcy, fr_i, zlight, zsoufer
      real(rk) :: chemo2, xfecolagg, zklight, consfe3, za1, plig, nitrfac, xcoagfe, zconsfe

      _LOOP_BEGIN_
         _GET_(self%id_fer, fer)
         _GET_(self%id_doc, doc)
         _GET_(self%id_poc, poc)
         _GET_(self%id_goc, goc)
         _GET_(self%id_cal, cal)
         _GET_(self%id_gsi, gsi)
         _GET_(self%id_no3, no3)
         _GET_(self%id_hi, hi)
         _GET_(self%id_oxy, oxy)
         _GET_(self%id_etot, etot)
         _GET_(self%id_xdiss, xdiss)
         _GET_SURFACE_(self%id_gphit, gphit)
         _GET_(self%id_tempis, tempis)
         _GET_(self%id_salinprac, salinprac)
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_zdust, zdust)
         _GET_(self%id_chemo2, chemo2)
         _GET_(self%id_consfe3_sum, consfe3)
         _GET_(self%id_nitrfac, nitrfac)

         ztkel = tempis + 273.15_rk
         zsal  = salinprac !(ji,jj,1) + ( 1.- tmask(ji,jj,1) ) * 35.
         zis    = 19.924_rk * zsal / ( 1000._rk- 1.005_rk * zsal )

         fekeq  = 10**( 17.27 - 1565.7 / ztkel )

         ! Liu and Millero (1999) only valid 5 - 50 degC
         ztkel1 = MAX( 5._rk , tempis ) + 273.16_rk
         fesol(1) = 10**(-13.486_rk - 0.1856_rk* zis**0.5 + 0.3073_rk*zis + 5254.0_rk/ztkel1)
         fesol(2) = 10**(2.517_rk - 0.8885_rk*zis**0.5 + 0.2139_rk * zis - 1320.0_rk/ztkel1 )
         fesol(3) = 10**(0.4511_rk - 0.3305_rk*zis**0.5 - 1996.0_rk/ztkel1 )
         fesol(4) = 10**(-0.2965_rk - 0.7881_rk*zis**0.5 - 4086.0_rk/ztkel1 )
         fesol(5) = 10**(4.4466_rk - 0.8505_rk*zis**0.5 - 7980.0_rk/ztkel1 )

         ! Total ligand concentration : Ligands can be chosen to be constant or variable
         ! Parameterization from Tagliabue and Voelker (2011)
         ! -------------------------------------------------

         xfecolagg = self%ligand * 1E9_rk + MAX(0._rk, chemo2 - oxy ) / 400.E-6_rk

         IF( ln_ligand ) THEN  ; !ztotlig = lgw * 1E9      ! Jorn: TODO
         ELSE
            IF( ln_ligvar ) THEN
               ztotlig =  0.09_rk * 0.667_rk * doc * 1E6_rk + xfecolagg
               ztotlig =  MIN( ztotlig, 10._rk )
            ELSE
               ztotlig = self%ligand * 1E9
            ENDIF
         ENDIF
          

         ! ------------------------------------------------------------
         !  from Aumont and Bopp (2006)
         ! This model is based on one ligand and Fe' 
         ! Chemistry is supposed to be fast enough to be at equilibrium
         ! ------------------------------------------------------------
         !consfe3 = 0._rk

         zTL1  = ztotlig
         zkeq            = fekeq
         zklight         = 4.77E-7_rk * etot * 0.5_rk / ( 10**(-6.3) )
         zconsfe =  consfe3 / ( 10**(-6.3) ) 
         zfesatur        = zTL1 * 1E-9_rk
         ztfe            = (1.0_rk + zklight) * fer

         ! Fe' is the root of a 2nd order polynom
         za1 =  1. + zfesatur * zkeq + zklight +  zconsfe - zkeq * fer


         _SET_DIAGNOSTIC_(self%id_zkeq_diag, zkeq)
         _SET_DIAGNOSTIC_(self%id_zklight_diag, zklight)
         _SET_DIAGNOSTIC_(self%id_zconsfe_diag, zconsfe)
         _SET_DIAGNOSTIC_(self%id_za1_diag, za1)
         _SET_DIAGNOSTIC_(self%id_ztfe_diag, ztfe)
          
         zFe3  = ( -1 * za1 + SQRT( za1**2 + 4._rk * ztfe * zkeq) ) / (2._rk * zkeq + rtrn ) ! Eq 65
         
         zFeL1 = MAX( 0._rk, fer - zFe3 )  ! Jorn: "complexed" iron (nmol/L)
         
         _SET_DIAGNOSTIC_(self%id_Fe3, zFe3)
         _SET_DIAGNOSTIC_(self%id_FeL1, zFeL1)
         _SET_DIAGNOSTIC_(self%id_zTL1, zTL1)

         plig =  MAX( 0._rk, ( zFeL1 / ( fer + rtrn ) ) )

         _SET_DIAGNOSTIC_(self%id_plig, plig)

         IF (ln_ligand) THEN
           zfecoll = 0.5 * zFeL1 * MAX(0._rk, ztotlig - xfecolagg) / ( ztotlig + rtrn ) 
         ELSE
            IF(ln_ligvar) THEN
               zfecoll = 0.5 * zFeL1 * MAX(0._rk, ztotlig - xfecolagg ) / ( ztotlig + rtrn )
            ELSE
               zfecoll = 0.5 * zFeL1 * MAX(0._rk, 0.09_rk * 0.667_rk * doc * 1E6_rk )/ ( ztotlig + rtrn )
            ENDIF
         ENDIF

         ! Scavenging rate of iron. This scavenging rate depends on the load of
         ! particles of sea water. 
         ! This parameterization assumes a simple second order kinetics
         ! (k[Particles][Fe]).
         ! Scavenging onto dust is also included as evidenced from the DUNE
         ! experiments.
         ! --------------------------------------------------------------------------------------

         zhplus  = max( rtrn, hi )
         fe3sol  = fesol(1) * ( zhplus**3 + fesol(2) * zhplus**2  &
         &         + fesol(3) * zhplus + fesol(4)     &
         &         + fesol(5) / zhplus )
         !
         !zfeequi = zFe3 * 1E-9
         ! precipitation of Fe3+, creation of nanoparticles
         precip = MAX( 0._rk, ( zFe3  - fe3sol ) ) * self%kfep * xstep * ( 1.0 - nitrfac )   ! Jorn: replaces Eq 62?
         precipno3 = 2.0_rk * 130.0_rk * no3 * nitrfac * xstep * zFe3
         zfeprecip = precip + precipno3
         !
         ztrc   = MAX( ( poc + goc + cal + gsi ) * 1.e6_rk , rtrn)
         
         zxlam  = MAX( 1.E-3_rk, (1._rk - EXP(-2 * oxy / 100.E-6_rk )))
         zlam1b = 3.e-5_rk + ( self%xlamdust * zdust + self%xlam1 * ztrc ) * zxlam
         zscave = zFe3 * zlam1b * xstep
         
         zlam1a   = ( 12.0_rk  * 0.3_rk * doc + 9.05_rk  * poc ) * xdiss    &
                 &    + ( 2.49_rk  * poc )     &
                 &    + ( 127.8_rk * 0.3_rk * doc + 725.7_rk * poc )

         zaggdfea = zlam1a * xstep * zfecoll  !Eq 61a - but Brownian POC term (a4) missing
         zlam1b   = ( 1.94_rk * xdiss + 1.37_rk ) * goc
         zaggdfeb = zlam1b * xstep * zfecoll  ! Eq 61b
         
         xcoagfe = zlam1a + zlam1b


         !  Compute the coagulation of colloidal iron. This parameterization 
         !  could be thought as an equivalent of colloidal pumping.
         !  It requires certainly some more work as it is very poorly constrained.
         !  ----------------------------------------------------------------
         !
         !
         _ADD_SOURCE_(self%id_fer, - zscave - zaggdfea - zaggdfeb - zfeprecip)
         _ADD_SOURCE_(self%id_sfe, zscave * self%scaveff * poc / ztrc   )
         _ADD_SOURCE_(self%id_bfe, zscave * self%scaveff * poc / ztrc   )

         _ADD_SOURCE_(self%id_sfe,  + zaggdfea)
         _ADD_SOURCE_(self%id_bfe,   + zaggdfeb)

         ! Precipitated iron is supposed to be permanently lost.
         ! Scavenged iron is supposed to be released back to seawater
         ! when POM is solubilized. This is highly uncertain as probably
         ! a significant part of it may be rescavenged back onto 
         ! the particles. An efficiency factor is applied that is read
         ! in the namelist. 
         ! See for instance Tagliabue et al. (2019).
         ! Aggregated FeL is considered as biogenic Fe as it 
         ! probably remains  complexed when the particle is solubilized.
         _SET_DIAGNOSTIC_(self%id_scav, 1.e9 * zscave/xstep)
         _SET_DIAGNOSTIC_(self%id_coll, 1.e9 * (zaggdfea + zaggdfeb)/xstep)

         _SET_DIAGNOSTIC_(self%id_zdust_diag, zdust)
         !
         !
         !  Define the bioavailable fraction of iron
         !  ----------------------------------------
         !biron = fer   ! Jorn: unused????

         !IF( ln_ligand ) THEN   ! Jorn: TODO, ln_ligand not implemented yet
         !   !
         !   zlam1a   = ( 0.369  * 0.3 * doc + 102.4  * poc ) * xdiss    &
         !         &    + ( 114.   * 0.3 * doc )
         !   !
         !   zlam1b   = 3.53E3 *   goc * xdiss
         !   zligco   = 0.5 * lgw
         !   zaggliga = zlam1a * xstep * zligco
         !   zaggligb = zlam1b * xstep * zligco
         !   _ADD_SOURCE_(self%id_lgw, - zaggliga - zaggligb)
         !   !zlcoll3d(ji,jj,jk)  = zaggliga + zaggligb
         !   !
         !   plig =  MAX( 0., ( ( zFeL1 * 1E-9 ) / ( fer +rtrn ) ) )
         !   !
         !ENDIF

         ! Additional non-conservative iron production controlled by light, e.g., photolysis, OA 2021-09-03 [from p4zsed, nitrogen fixation section]
         _GET_(self%id_etot_ndcy, etot_ndcy)
         _GET_SURFACE_(self%id_fr_i, fr_i)
         zlight  =  ( 1._rk- EXP( -etot_ndcy / self%light ) ) * ( 1. - fr_i )  ! Jorn: light limitation of diazotrophs (Eq 58b), reused here for iron photolysis
         zsoufer = zlight * 2E-11_rk / ( 2E-11_rk + fer )                         ! Jorn: limitation factor that is 1 when ambient iron is zero, and 0 in iron-replete environments
         _ADD_SOURCE_(self%id_fer, + 0.002_rk * 4E-10_rk * zsoufer / rday)        ! Jorn : dropped multiplication with rfact2 [time step in seconds]

      _LOOP_END_
   end subroutine

end module
