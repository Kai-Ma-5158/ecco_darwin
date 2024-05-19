#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SURFFORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SURFFORCING(
     O           gDIC, gALK, gO2,
     O           gO2TH, gN2O, gN2OTH, gAR,!makai
     I           bi,bj,imin,imax,jmin,jmax,
     I           myTime,myIter,myThid)

C !DESCRIPTION:
C Update tendency terms for alkalinity, oxygen and DIC from air-sea
C exchanges

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#endif
#ifdef ALLOW_DARWIN
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_EXF_FIELDS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_FIELDS.h"
#include "DARWIN_FLUX.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  bi,bj                :: tile indices
      INTEGER iMin,iMax,jMin,jMax,bi,bj
      INTEGER myIter, myThid
      _RL myTime

C !OUTPUT PARAMETERS: ==================================================
C  gDIC :: DIC tendency due to air-sea exchange
C  gALK :: ALK tendency due to air-sea exchange
C  gO2  :: O2 tendency due to air-sea exchange
      _RL  gDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  gALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  gO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C !makai 
      _RL gO2TH(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL gN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL gN2OTH(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL gAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
C-!makai

#ifdef ALLOW_DARWIN
#ifdef DARWIN_ALLOW_CARBON

C !LOCAL VARIABLES: ====================================================
      INTEGER i,j,k,ks
      _RL co3dummy
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
      _RL SchmidtNoDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch0(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C !makai
C Adding  Sol. relationship coeff. for N2O andArgon
      _RL KwexchAr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL KwexchN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
C !makai

      _RL pisvel0(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C local variables for carbon chem
      _RL surfdic
      _RL surfalk
      _RL surfphos
      _RL surfsi
      _RL surftemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL SchmidtNoO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C makai 2024/03/08
C Adding here variables for the Schmidt Number
C of N2O and Argon
      _RL SchmidtNoAr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL SchmidtNoN2O(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
C makai
       _RL pN2Oatm !makai
C Setting atmospheric pressure of N2O!makai
C units are in ppb or nanoatm!makai
C example 320 ppb = 320E-9!makai
C         pN2Oatm = 320E-9!makai
C !makai


      _RL O2sat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ARsat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL N2Osat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)!makai
      _RL aTT
      _RL aTK
      _RL aTS
      _RL aTS2
      _RL aTS3
      _RL aTS4
      _RL aTS5
      _RL o2s
      _RL ttemp
      _RL stemp
      _RL oCnew
CEOP
C-makai
C NEW VARIABLES!makai
C FOR ARGON!makai
      _RL rhosw!makai
      _RL sumAr!makai
      _RL Arsol!makai
      _RL TempS !makai    
C N2O Variables!makai
      _RL tk!makai
      _RL soln2o!makai

      _RL a1,a2,a3,b1,b2,b3!makai
      _RL part1, part2, part3!makai
C      _RL sal!makai
C      _RL  pN2Oatm!makai

C coeffs for N2O !makai
         a1 = -62.7062!makai
         a2 =  97.3066!makai
         a3 =  24.1406!makai
         b1 = -0.05842!makai
         b2 =  0.033193!makai
         b3 = -0.0051313!makai
     
         pN2Oatm = 320E-9 !makai

C  Density of seawater (kg/m3)!makai
       rhosw = 1024.5!makai
CEOP
      ks = 1

      DO j=jmin,jmax
       DO i=imin,imax
         gDIC(i,j) = 0.0 _d 0
         gALK(i,j) = 0.0 _d 0
         gO2(i,j) = 0.0 _d 0
C-makai
C Four new tracers !makai
         gO2TH(i,j) = 0.0 _d 0!makai
         gN2O(i,j) = 0.0 _d 0 !makai  
         gN2OTH(i,j) = 0.0 _d 0!makai
         gAR(i,j) = 0.0 _d 0!makai
       ENDDO
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Compute AtmosP and Kwexch0 which are used for flux of CO2 and O2
      DO j=jmin,jmax
       DO i=imin,imax
        IF (maskC(i,j,ks,bi,bj).NE.0. _d 0) THEN

#ifdef USE_PLOAD
C Convert anomalous pressure pLoad (in Pa) from atmospheric model
C to total pressure (in Atm)
C Note: it is assumed the reference atmospheric pressure is 1Atm=1013mb
C       rather than the actual ref. pressure from Atm. model so that on
C       average AtmosP is about 1 Atm.
         AtmosP(i,j,bi,bj)= 1. _d 0 + pLoad(i,j,bi,bj)/Pa2Atm
#endif

C Pre-compute part of exchange coefficient: pisvel*(1-fice)
C Schmidt number is accounted for later
         pisvel0(i,j) = 0.337 _d 0 * windSpeed(i,j,bi,bj)**2/3.6 _d 5
         Kwexch0(i,j) = pisvel0(i,j) *
     &                   (1. _d 0 - iceFrac(i,j,bi,bj))

        ENDIF
       ENDDO
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c flux of alkalinity

#ifdef ALLOW_OLD_VIRTUALFLUX
      DO j=jmin,jmax
       DO i=imin,imax
        IF (maskC(i,j,ks,bi,bj).NE.0. _d 0) THEN
c calculate virtual flux
c EminusPforV = dS/dt*(1/Sglob)
C NOTE: Be very careful with signs here!
C Positive EminusPforV => loss of water to atmos and increase
C in salinity. Thus, also increase in other surface tracers
C (i.e. positive virtual flux into surface layer)
         VFluxAlk(i,j,bi,bj) = gsm_ALK*surfaceForcingS(i,j,bi,bj)/gsm_S
         gALK(i,j) =
     &     recip_drF(ks) * VFluxAlk(i,j,bi,bj)
        ELSE
         VFluxAlk(i,j,bi,bj) = 0 _d 0
        ENDIF
       ENDDO
      ENDDO
#endif /* ALLOW_OLD_VIRTUALFLUX */

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c flux of O2

C calculate SCHMIDT NO. for O2
      DO j=jmin,jmax
       DO i=imin,imax
        IF (maskC(i,j,ks,bi,bj).NE.0.) THEN
         ttemp = theta(i,j,ks,bi,bj)
         stemp = salt(i,j,ks,bi,bj)

         SchmidtNoO2(i,j) =
     &       sox1
     &     + sox2 * ttemp
     &     + sox3 * ttemp*ttemp
     &     + sox4 * ttemp*ttemp*ttemp

C-makai
C Computing now the Schmidt Number for Argon!makai
           SchmidtNoAr(i,j) =                          !makai
     &            sar1                                 !makai
     &          + sar2 * ttemp                         !makai
     &          + sar3 * ttemp*ttemp                   !makai
     &          + sar4 * ttemp*ttemp*ttemp             !makai
C        write(*,*)'SchmidtNoAr = ',  SchmidtNoAr(i,j) !makai   

C Computing now the Schmidt number for N2O!makai

        SchmidtNoN2O(i,j) =                             !makai
     &            sn2o1                                 !makai
     &          + sn2o2 * ttemp                         !makai
     &          + sn2o3 * ttemp*ttemp                   !makai
     &          + sn2o4 * ttemp*ttemp*ttemp             !makai
C       write(*,*)'SchmidtNoN2O(i,j) =', SchmidtNoN2O(i,j)!makai
C-makai

C Determine surface flux of O2
C exchange coeff accounting for ice cover and Schmidt no.
C Kwexch0= pisvel*(1-fice): previously computed in dic_surfforcing.F

         Kwexch(i,j) = Kwexch0(i,j)
     &                / sqrt(SchmidtNoO2(i,j)/660.0 _d 0)

C-makai
C convert T to scaled temperature for Argon calculation    !makai
           tempS = log((298.15 - ttemp)/(273.15 + ttemp)); !makai
          KwexchAr(i,j) = Kwexch0(i,j)                     !makai
     &                / sqrt(SchmidtNoAr(i,j)/660.0 _d 0)  !makai
C            write(*,*)'KwexchAr =', KwexchAr(i,j)         !makai

C Now same for N2O                                         !makai

          KwexchN2O(i,j) = Kwexch0(i,j)                    !makai
     &                / sqrt(SchmidtNoN2O(i,j)/660.0 _d 0) !makai
C           write(*,*)'KwexchN2O =', KwexchN2O(i,j)        !makai
C determine saturation Ar                                  !makai
C using Hamme and Emerson                                  !makai
          sumAr=  A0ar+(A1ar*tempS)+(A2ar*tempS**2)+       !makai
     &   (A3ar*(tempS**3))+stemp*(B0ar+(B1ar*tempS)+       !makai
     &   (B2ar*(tempS**2)))                                !makai

C  Solubilty factor computed in micromol/Kg                !makai
C  and conversion into mmol/m3                             !makai

           Arsol=EXP(sumAr)*(rhosw/1e3)      !makai
C conversion into mol/m3                     !makai
C           Arsol=Arsol/1e3                  !makai
          Arsat(i,j)=Arsol                   !makai
C          write(*,*)'Arsat =', Arsat(i,j)   !makai
C-makai

C determine saturation O2
C using Garcia and Gordon (1992), L&O (mistake in original ?)
         aTT  = 298.15 _d 0 -ttemp
         aTK  = 273.15 _d 0 +ttemp
         aTS  = log(aTT/aTK)
         aTS2 = aTS*aTS
         aTS3 = aTS2*aTS
         aTS4 = aTS3*aTS
         aTS5 = aTS4*aTS

         oCnew = oA0 + oA1*aTS + oA2*aTS2 + oA3*aTS3 +
     &           oA4*aTS4 + oA5*aTS5 +
     &           (oB0 + oB1*aTS + oB2*aTS2 + oB3*aTS3)*stemp +
     &           oC0*stemp*stemp

         o2s = EXP(oCnew)

c molar volume of O2: O2mol2L = 22.391 L mol-1 = 0.022391 L mmol-1
c o2s in ml/l = l/m3
c O2sat = o2s / O2mmol2L  (in mmol/m3)
c Convert from ml/l to mmol/m^3
         O2sat(i,j) = o2s/22.3916 _d -3

C Determine flux, inc. correction for local atmos surface pressure
         FluxO2(i,j,bi,bj) = Kwexch(i,j)*
     &                     (AtmosP(i,j,bi,bj)*O2sat(i,j)
     &                      - pTracer(i,j,ks,bi,bj,iO2))

C-makai
C For O2 Thermal Only (O2TH)                                  !makai
       FluxO2TH(i,j,bi,bj) = Kwexch(i,j)*                     !makai
     &                     (AtmosP(i,j,bi,bj)*O2sat(i,j)      !makai
     &                      - pTracer(i,j,ks,bi,bj,iO2TH))    !makai     
        
C Determine flux, inc. correction for local atmos surface pressure !makai
              FluxAR(i,j,bi,bj) = KwexchAr(i,j)*             !makai
     &                     (AtmosP(i,j,bi,bj)*Arsat(i,j)     !makai
     &                      - pTracer(i,j,ks,bi,bj,iAR))     !makai


C Now flux of N2O and N2O-TH                                  !makai
C Convert celsius into kelvin                                 !makai
             tk = ttemp+273.15                                !makai
C    Coefficients of solubilty factors of N2O                 !makai
            soln2o=(exp(a1+a2*(100/tk)+a3*log(tk/100)         !makai
     &     +stemp*(b1+b2*(tk/100)+b3*((tk/100)**2))))         !makai  
C convert into mol/m3/atm                                     !makai
             soln2o = soln2o*1E3                              !makai
c Determine flux, inc. correction for local atmos surface pressure
cQQ PTR_O2?
C The flux is calculated in mmol N2O
C - N2O TOTAL
              FluxN2O(i,j,bi,bj) = KwexchN2O(i,j)*            !makai
     &                     ((pN2Oatm*soln2o*1E3)              !makai
     &                      -  pTracer(i,j,ks,bi,bj,iN2O))    !makai
C - N2O THERMAL ONLY
              FluxN2OTH(i,j,bi,bj) = KwexchN2O(i,j)*           !makai
     &                     ((pN2Oatm*soln2o*1E3)               !makai
     &                      -  pTracer(i,j,ks,bi,bj,iN2OTH))   !makai
C-makai

C Saving tendencies here
         gO2(i,j) =
     &      recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*FluxO2(i,j,bi,bj)
C-makai
         gO2TH(i,j) =                                              !makai
     &      recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*FluxO2TH(i,j,bi,bj)!makai
         gAR(i,j) =                                                    !makai
     &      recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*FluxAR(i,j,bi,bj) !makai
     
         gN2O(i,j) =                                                    !makai
     &      recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*FluxN2O(i,j,bi,bj) !makai

         gN2OTH(i,j) =                                                 !makai
     &      recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*FluxN2OTH(i,j,bi,bj)!makai
C-makai
        ENDIF
       ENDDO
      ENDDO

C ======================================================================
      DO k=1,Nr
C ======================================================================

C determine inorganic carbon chem coefficients
      DO j=jmin,jmax
       DO i=imin,imax
c put bounds on tracers so pH solver doesn't blow up
        surfsalt(i,j) = MAX(surfSaltMin,
     &                  MIN(surfSaltMax, salt(i,j,k,bi,bj)))
        surftemp(i,j) = MAX(surfTempMin,
     &                  MIN(surfTempMax, theta(i,j,k,bi,bj)))
       ENDDO
      ENDDO

      CALL DARWIN_CARBON_COEFFS(
     I                   surftemp,surfsalt,
     I                   bi,bj,iMin,iMax,jMin,jMax,k,myThid)
C====================================================================

c pCO2 solver...
      DO j=jmin,jmax
       DO i=imin,imax
        IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) THEN
c put bounds on tracers so pH solver doesn't blow up
         surfdic  = ptr2mol * MAX(surfDICMin,
     &              MIN(surfDICMax, Ptracer(i,j,k,bi,bj,iDIC)))
     &              * maskC(i,j,k,bi,bj)
         surfalk  = ptr2mol * MAX(surfALKMin,
     &              MIN(surfALKMax, Ptracer(i,j,k,bi,bj,iALK)))
     &              * maskC(i,j,k,bi,bj)
         surfphos = ptr2mol * MAX(surfPO4Min,
     &              MIN(surfPO4Max, Ptracer(i,j,k,bi,bj,iPO4)))
     &              * maskC(i,j,k,bi,bj)
         surfsi   = ptr2mol * MAX(surfSiMinINit,
     &              MIN(surfSiMax, Ptracer(i,j,k,bi,bj,iSiO2)))
     &              * maskC(i,j,k,bi,bj)
         CALL DARWIN_CALC_PCO2_APPROX(
     I        surftemp(i,j), surfsalt(i,j),
     I        surfdic, surfphos, surfsi, surfalk,
     I        ak1(i,j,bi,bj), ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj), ak2p(i,j,bi,bj), ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj), akb(i,j,bi,bj), akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj), akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj),  fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj), st(i,j,bi,bj), ft(i,j,bi,bj),
     U        pH(i,j,k,bi,bj), pCO2(i,j,k,bi,bj), co3dummy,
     I        i,j,k,bi,bj,myIter,myThid )
        ELSE
         pCO2(i,j,k,bi,bj)=0. _d 0
        ENDIF
       ENDDO
      ENDDO

      IF (k .EQ. ks) THEN

      DO j=jmin,jmax
       DO i=imin,imax
        IF ( maskC(i,j,ks,bi,bj).NE.0. _d 0 ) THEN
         ttemp = theta(i,j,ks,bi,bj)
C calculate SCHMIDT NO. for CO2
         SchmidtNoDIC(i,j) =
     &      sca1 +
     &      sca2*ttemp +
     &      sca3*ttemp*ttemp +
     &      sca4*ttemp*ttemp*ttemp
c make sure Schmidt number is not negative (will happen if temp>39C)
         SchmidtNoDIC(i,j) = MAX(1.0 _d -2, SchmidtNoDIC(i,j))

C Determine surface flux (FDIC)
C first correct apCO2 for surface atmos pressure
         apCO2Sat(i,j,bi,bj) = AtmosP(i,j,bi,bj)*atmospCO2(i,j,bi,bj)

C then account for Schmidt number
         Kwexch(i,j) = Kwexch0(i,j)
     &                 / sqrt(SchmidtNoDIC(i,j)/660.0 _d 0)

C compute fugacity of CO2     
         fugCO2(i,j,bi,bj) = 
     &    pCO2(i,j,ks,bi,bj)*fugf(i,j,bi,bj)
     
c Flux = kw*rho*(ff*pCO2atm-k0*FugFac*pCO2ocean)
         FluxCO2(i,j,bi,bj) =
     &          Kwexch(i,j)*(
     &            ff(i,j,bi,bj)*apCO2Sat(i,j,bi,bj) -
     &            pCO2(i,j,ks,bi,bj)*fugf(i,j,bi,bj)*ak0(i,j,bi,bj) )
        ELSE
         FluxCO2(i,j,bi,bj) = 0. _d 0
        ENDIF
C convert flux (mol kg-1 m s-1) to (mmol m-2 s-1)
        FluxCO2(i,j,bi,bj) = FluxCO2(i,j,bi,bj)/m3perkg/ptr2mol

#ifdef ALLOW_OLD_VIRTUALFLUX
        IF (maskC(i,j,ks,bi,bj).NE.0. _d 0) THEN
c calculate virtual flux
c EminusPforV = dS/dt*(1/Sglob)
C NOTE: Be very careful with signs here!
C Positive EminusPforV => loss of water to atmos and increase
C in salinity. Thus, also increase in other surface tracers
C (i.e. positive virtual flux into surface layer)
C ...so here, VirtualFLux = dC/dt!
         VFluxCO2(i,j,bi,bj)=gsm_DIC*surfaceForcingS(i,j,bi,bj)/gsm_S
c OR
c let virtual flux be zero
c        VirtualFlux(i,j)=0.d0
c
        ELSE
         VFluxCO2(i,j,bi,bj)=0. _d 0
        ENDIF
#endif /* ALLOW_OLD_VIRTUALFLUX */
       ENDDO
      ENDDO

C update tendency
      DO j=jmin,jmax
       DO i=imin,imax
        IF ( maskC(i,j,ks,bi,bj).NE.0. _d 0 ) THEN
         gDIC(i,j) =
     &     recip_drF(ks)*recip_hFacC(i,j,ks,bi,bj)*
     &              ( FluxCO2(i,j,bi,bj)
#ifdef ALLOW_OLD_VIRTUALFLUX
     &              + VFluxCO2(i,j,bi,bj)
#endif
     &              )
        ENDIF
       ENDDO
      ENDDO

      CALL DIAGNOSTICS_FILL(fugf,    'fugfCO2 ',0,1,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(fugCO2,  'fCO2    ',0,1,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(FluxO2,  'fluxO2  ',0,1,1,bi,bj,myThid)
#ifdef ALLOW_OLD_VIRTUALFLUX
      CALL DIAGNOSTICS_FILL(VFluxCO2,'VfluxCO2',0,1,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(VFluxAlk,'VfluxAlk',0,1,1,bi,bj,myThid)
#endif

C     k is at surface
      ENDIF

C === k ================================================================
      ENDDO
C ======================================================================

      CALL DIAGNOSTICS_FILL(atmospCO2,'apCO2   ',0,1,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(apCO2Sat,'apCO2sat',0,1,1,bi,bj,myThid)

#endif /* DARWIN_ALLOW_CARBON */
#endif /* ALLOW_DARWIN */

      RETURN
      END
