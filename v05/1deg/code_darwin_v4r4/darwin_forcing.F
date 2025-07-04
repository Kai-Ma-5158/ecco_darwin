C TODO
C - have to call fechem after model?

#include "DARWIN_OPTIONS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_OPTIONS.h"
#endif

CBOP
C !ROUTINE: DARWIN_FORCING
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_FORCING( Ptrdummy,
     I                        bi, bj, iMin, iMax, jMin, jMax,
     I                        myIter, myTime, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE.h"
#endif
#ifdef ALLOW_RADTRANS
#include "RADTRANS_SIZE.h"
#endif
#ifdef ALLOW_DARWIN
#include "GCHEM.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_DIAGS.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"
#include "DARWIN_FIELDS.h"
#include "DARWIN_EXF_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  Ptrdummy             :: dummy for darwin2 compatibility
C  myThid               :: thread number
      _RL Ptrdummy(*)
      _RL myTime
      INTEGER iMin, iMax, jMin, jMax, bi, bj, myIter, myThid
CEOP

#ifdef ALLOW_DARWIN

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  k                    :: vertical level
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      INTEGER i,j,k,kdn,iTr,l,isub
      CHARACTER*8 diagname
      _RL dTsub(Nr)
C-Makai
C-array of values of depth of tracers
      _RL depthT(Nr)
C-Makai

      _RL midTime, subTime
      _RL gPtr(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nDarwin)
      _RL gPtrTemp(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nDarwin)
      _RL PAR(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nlam)
      _RL diags(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, darwin_nDiag)
      _RL gDIC(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gALK(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gDICsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gALKsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNO3surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNO2surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gNH4surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)  
      _RL gPO4surfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gFeTsurfForc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL gSiO2surfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gO2surfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gO2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
C-Makai
      _RL gO2THsurfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gO2TH(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
C
      _RL gN2OsurfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gN2O(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
C

      _RL gN2OTHsurfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gN2OTH(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
C

      _RL gARsurfForc(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL gAR(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)

C-Makai
      _RL surfFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL freeFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe1(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL diffFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
#ifndef DARWIN_ALLOW_CHLQUOTA
      _RL tmp3d(1:sNx,1:sNy,Nr)
#endif
      _RL scv,scav_poc
      _RL flx, POCl
      _RL ptr(nDarwin), gtr(nDarwin), PARl(nlam)

C PAR in the first layer - surface
      _RL PARlsrf(nlam)!makai
C and depth local var
      _RL depthl!makai

      _RL chlout(nPhoto)
      _RL diagsl(darwin_nDiag)
      _RL photoTempFunc(nplank)
      _RL hetTempFunc(nplank)
      _RL grazTempFunc(nplank)
      _RL reminTempFunc
      _RL mortTempFunc
      _RL mort2TempFunc
      _RL uptakeTempFunc
      _RL tmp
      _RL sedFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cNH4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)      
      _RL cPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cSi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sNH4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sSi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sFe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL consDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL consDIC_PIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL respirDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL reminDIC_DOC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL reminDIC_POC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL disscDIC_PIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL cO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

      DO k=1,Nr
       dTsub(k) = PTRACERS_dTLev(k)/nsubtime
       depthT(k) = drC(k)!makai
      ENDDO

C     time at middle of sub-timestep
      midTime = myTime - deltaTclock + .5*deltaTclock/nsubtime
C     time at end of sub-timestep
      subTime = myTime - deltaTclock + deltaTclock/nsubtime
      DO isub=1,nsubtime

C === reset tendencies =================================================
      DO itr=1,nDarwin
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        gPtr(i,j,k,iTr) = 0.0 _d 0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        diffFe1(i,j,k) = 0.0 _d 0
        diffFe2(i,j,k) = 0.0 _d 0
        diffFe(i,j,k) = 0.0 _d 0
      ENDDO
      ENDDO
      ENDDO
C === light ============================================================
C     Initialize Chl from balanced-growth Chl:C if requested
C     and check Chl:C bounds.
C     Note: myIter has already been incremented
c      IF (myIter-1 .EQ. darwin_chlIter0) THEN
c        CALL DARWIN_INIT_CHL(bi, bj, subTime, myIter, myThid)
c      ENDIF

#ifdef ALLOW_SEAICE
      IF (DARWIN_useSEAICE) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         iceFrac(i,j,bi,bj) = AREA(i,j,bi,bj)
        ENDDO
       ENDDO
      ENDIF
#endif

      CALL TIMER_START('DARWIN_LIGHT [DARWIN_FORCING]',myThid)
#ifdef ALLOW_RADTRANS
      CALL DARWIN_LIGHT_RADTRANS(PAR,subTime,bi,bj,iMin,iMax,jMin,jMax,
     &                        subTime,myIter,myThid)
#else
      CALL DARWIN_LIGHT(PAR, midTime, bi, bj, iMin, iMax, jMin, jMax, 
     &               subTime, myIter, myThid)
#endif
      CALL TIMER_STOP ('DARWIN_LIGHT [DARWIN_FORCING]',myThid)

C === dic ==============================================================
#ifdef DARWIN_ALLOW_CARBON
C carbon air-sea interaction
      CALL TIMER_START('DARWIN_SURFFORCING [DARWIN_FORCING]',myThid)
      CALL DARWIN_SURFFORCING(
     O                    gDIC, gALK, gO2,gO2TH,gN2O,gN2OTH,gAR,!makai
     I                    bi,bj,iMin,iMax,jMin,jMax,
     I                    subTime,myIter,myThid)
      CALL TIMER_STOP ('DARWIN_SURFFORCING [DARWIN_FORCING]',myThid)
      DO j=jMin,jMax
      DO i=iMin,iMax
        gPtr(i,j,1,iDIC) = gPtr(i,j,1,iDIC) + gDIC(i,j)
        gPtr(i,j,1,iALK) = gPtr(i,j,1,iALK) + gALK(i,j)
        gPtr(i,j,1,iO2)  = gPtr(i,j,1,iO2)  + gO2(i,j)
C        gDICsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iDIC) 
C     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)	
      ENDDO
      ENDDO
#endif

C surface fluxes from E/P/runoff
      DO j=jMin,jMax
      DO i=iMin,iMax
        gDICsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iDIC)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gALKsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iALK)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNO3surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNO3)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gNH4surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iNH4)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gPO4surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iPO4)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gFeTsurfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iFeT)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gSiO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iSiO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
        gO2surfForc(i,j) = surfaceForcingPTr(i,j,bi,bj,iO2)
     &   * recip_drF(1)*recip_hFacC(i,j,1,bi,bj)
C-makai
C  ADDING FOUR NEW TRACERS
        gPtr(i,j,1,iO2TH)  = gPtr(i,j,1,iO2TH)  + gO2TH(i,j)!makai
        gPtr(i,j,1,iN2O)  = gPtr(i,j,1,iN2O)  + gN2O(i,j)!makai
        gPtr(i,j,1,iN2OTH)  = gPtr(i,j,1,iN2OTH)  + gN2OTH(i,j)!makai
        gPtr(i,j,1,iAR)  = gPtr(i,j,1,iAR)  + gAR(i,j)!makai
C-Makai

      ENDDO
      ENDDO

C === iron =============================================================
#ifdef DARWIN_ALLOW_CONS
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         DARWIN_minFeLoss(i,j,k,bi,bj) = Ptracer(i,j,k,bi,bj,iFeT)
        ENDDO
       ENDDO
      ENDDO
#endif
      CALL TIMER_START('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
      CALL DARWIN_FE_CHEM(
     U                 Ptracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 diffFe1(1-OLx,1-OLy,1),
     O                 freeFe(1-OLx, 1-OLy, 1),
     I                 bi, bj, iMin, iMax, jMin, jMax, myThid)
      CALL TIMER_STOP ('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
#ifdef DARWIN_ALLOW_CONS
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         DARWIN_minFeLoss(i,j,k,bi,bj) = DARWIN_minFeLoss(i,j,k,bi,bj)
     &                                 - Ptracer(i,j,k,bi,bj,iFeT)
        ENDDO
       ENDDO
      ENDDO
#endif

C     iron dust input
      CALL TIMER_START('DARWIN_FE_DUST [DARWIN_FORCING]',myThid)
      DO j=jMin,jMax
      DO i=iMin,iMax
      IF (hFacC(i,j,1,bi,bj) .NE. 0.) THEN
        gPtr(i,j,1,iFeT) =
     &  gPtr(i,j,1,iFeT) + alpfe * inputFe(i,j,bi,bj) /
     &                           drF(1) / hFacC(i,j,1,bi,bj)
      ENDIF
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_DUST [DARWIN_FORCING]',myThid)

C     scavenging
      CALL TIMER_START('DARWIN_FE_SCAV [DARWIN_FORCING]',myThid)
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
#ifdef DARWIN_PART_SCAV_POP
        scav_poc = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOP))/scav_R_POPPOC
        scv = scav_rat*scav_inter*(scav_poc**scav_exp)
#elif defined(DARWIN_PART_SCAV)
        scav_poc = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOC))
        scv = scav_rat*scav_inter*(scav_poc**scav_exp)
#else
        scv = scav
#endif
        gPtr(i,j,k,iFeT) = gPtr(i,j,k,iFeT) -
     &                     scv*MAX(0 _d 0, freefe(i,j,k))
#ifdef DARWIN_ALLOW_CONS
        DARWIN_partScav(i,j,k,bi,bj) = scv*MAX(0 _d 0, freefe(i,j,k))
#endif
      ENDDO
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_SCAV [DARWIN_FORCING]',myThid)

C     iron sediment source (in bottom grid cell above kMaxFeSed)
      CALL TIMER_START('DARWIN_FE_SED [DARWIN_FORCING]',myThid)
#ifdef DARWIN_ALLOW_CONS
      DO j=jMin,jMax
       DO i=iMin,iMax
        ironSedFlux(i,j,bi,bj) = 0 _d 0
       ENDDO
      ENDDO
#endif
      DO k = kMinFeSed, kMaxFeSed
       kdn = MIN(Nr, k+1)
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF (hFacC(i,j,k,bi,bj) .GT. 0. .AND.
     &       (hFacC(i,j,kdn,bi,bj) .EQ. 0. .OR. k .EQ. Nr)) THEN
#ifdef DARWIN_IRON_SED_SOURCE_VARIABLE
# ifdef DARWIN_IRON_SED_SOURCE_POP
          flx = fesedflux_pcm*wp_sink*R_CP_fesed*
     &            MAX(0 _d 0, Ptracer(i,j,k-1,bi,bj,iPOP))
# else
          POCl = MAX(0 _d 0, Ptracer(i,j,k,bi,bj,iPOC))
          flx = fesedflux_pcm*wc_sink*POCl - fesedflux_min
          flx = MAX(0. _d 0, flx)
# endif
#else
          flx = fesedflux
#endif
	  sedFe(i,j,k) = flx/(drF(k)*hFacC(i,j,k,bi,bj))
          gPtr(i,j,k,iFeT) = gPtr(i,j,k,iFeT) +
     &                             flx/(drF(k)*hFacC(i,j,k,bi,bj))
#ifdef DARWIN_ALLOW_CONS
          ironSedFlux(i,j,bi,bj) = flx
#endif
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_FE_SED [DARWIN_FORCING]',myThid)

C === plankton =========================================================
      CALL TIMER_START('DARWIN_PLANKTON [DARWIN_FORCING]',myThid)
      do k=1,Nr
C saving temp depth local value
      depthl=RC(k)!makai
      do j=jMin,jMax
      do i=iMin,iMax
      if (hFacC(i,j,k,bi,bj) .GT. 0.) then
        CALL DARWIN_TEMPFUNC(Theta(i,j,k,bi,bj),
     &         photoTempFunc, hetTempFunc, grazTempFunc,
     &         reminTempFunc, mortTempFunc, mort2TempFunc,
     &         uptakeTempFunc, myThid)
        DO iTr=1,nDarwin
         ptr(iTr) = Ptracer(i, j, k, bi, bj, iTr)
         gtr(iTr) = gPtr(i, j, k, iTr)
        ENDDO
        DO l=1,nlam
         PARl(l) = MAX(0 _d 0, PAR(i, j, k, l))
C-ADDing here the surface value only!makai
         PARlsrf(l)  = MAX(0 _d 0, PAR(i, j, 1, l))!makai
        ENDDO
        CALL DARWIN_PLANKTON(
     I                 ptr,
     U                 gtr,
     O                 chlout, diagsl,
     I                 PARl, photoTempFunc,
     I                 PARlsrf,depthl,!makai
     I                 hetTempFunc,
     I                 grazTempFunc,
     I                 reminTempFunc,
     I                 mortTempFunc, mort2TempFunc,
     I                 uptakeTempFunc,
#ifdef DARWIN_DEBUG
     I                 myxgloballo+(bi-1)*sNx+i-1,
     I                 myygloballo+(bj-1)*sNy+j-1,
     I                 k, dTsub(k),
#endif
     I                 subTime, myIter, myThid )
        DO l=1,darwin_nDiag
         diags(i, j, k, l) = diagsl(l)
        ENDDO
#ifdef DARWIN_ALLOW_CONS
        DARWIN_Nfix(i,j,k,bi,bj) = diagsl(iNfix)
        DARWIN_Ndenit(i,j,k,bi,bj) = diagsl(iDenitN)
#endif
#ifndef DARWIN_ALLOW_CHLQUOTA
#ifdef ALLOW_RADTRANS
        DO l=1,nPhoto
         chlPrev(i, j, k, bi, bj, l) = chlout(l)
        ENDDO
#else
        tmp = 0 _d 0
        DO l=1,nPhoto
         tmp = tmp + chlout(l)
        ENDDO
        chlPrev(i, j, k, bi, bj) = tmp
#endif
#endif
        DO iTr=1,nDarwin
         gPtr(i, j, k, iTr) = gtr(iTr)
        ENDDO
      endif
      enddo
      enddo
      enddo
      CALL TIMER_STOP('DARWIN_PLANKTON [DARWIN_FORCING]',myThid)

C === pre-sinking diags ================================================
#ifdef ALLOW_DIAGNOSTICS
      IF (useDIAGNOSTICS) THEN
       CALL TIMER_START('DIAGS_FILL [DARWIN_FORCING]',myThid)
       DO iTr=1,nDarwin
        diagname = 'gECO'//PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL(gPtr(1-OLx,1-OLy,1,iTr), diagname,
     &          0,Nr,2,bi,bj,myThid)
       ENDDO
       CALL TIMER_STOP ('DIAGS_FILL [DARWIN_FORCING]',myThid)
      ENDIF
#endif

C === sinking ==========================================================
      CALL TIMER_START('DARWIN_SINKING [DARWIN_FORCING]',myThid)
      CALL DARWIN_SINKING( Ptracer,gPtr,bi,bj,iMin,iMax,jMin,jMax,
     &                  subTime,myIter,myThid )
      CALL TIMER_STOP ('DARWIN_SINKING [DARWIN_FORCING]',myThid)

C === apply tendencies to tracers ======================================
      CALL TIMER_START('DARWIN_STEP [DARWIN_FORCING]',myThid)
      DO iTr=1,nDarwin
      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        pTracer(i,j,k,bi,bj,iTr)=pTracer(i,j,k,bi,bj,iTr)
     &                      +dTsub(k)*gPtr(i,j,k,iTr)*maskInC(i,j,bi,bj)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL TIMER_STOP ('DARWIN_STEP [DARWIN_FORCING]',myThid)

C === iron =============================================================
C     re-apply free iron limit to FeT
#ifdef DARWIN_ALLOW_CONS
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         DARWIN_minFeLoss(i,j,k,bi,bj) = DARWIN_minFeLoss(i,j,k,bi,bj)
     &                                 + Ptracer(i,j,k,bi,bj,iFeT)
        ENDDO
       ENDDO
      ENDDO
#endif
      CALL TIMER_START('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
      CALL DARWIN_FE_CHEM(
     U                 pTracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 diffFe2(1-OLx,1-OLy,1),
     O                 freeFe(1-OLx,1-OLy,1),
     I                 bi, bj, iMin, iMax, jMin, jMax, myThid)
      CALL TIMER_STOP ('DARWIN_FE_CHEM [DARWIN_FORCING]',myThid)
#ifdef DARWIN_ALLOW_CONS
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         DARWIN_minFeLoss(i,j,k,bi,bj) = DARWIN_minFeLoss(i,j,k,bi,bj)
     &                                 - Ptracer(i,j,k,bi,bj,iFeT)
        ENDDO
       ENDDO
      ENDDO
#endif

      DO k=1,Nr
      DO j=jMin,jMax
      DO i=iMin,iMax
        diffFe(i,j,k) = hFacC(i,j,k,bi,bj) 
     &   * ((diffFe1(i,j,k) + diffFe2(i,j,k)) / deltaTClock)
      ENDDO
      ENDDO
      ENDDO
C === diagnostics ======================================================
#ifdef ALLOW_DIAGNOSTICS
      IF (useDIAGNOSTICS) THEN
       CALL TIMER_START('DIAGS_FILL [DARWIN_FORCING]',myThid)
       DO l = 1, nlam
        WRITE(diagname, '(A,I3.3)') 'PAR', l
        CALL DIAGNOSTICS_FILL(PAR(1-OLx,1-OLy,1,l),diagname,0,Nr,2,
     &        bi,bj,myThid)
       ENDDO
       IF (DIAGNOSTICS_IS_ON('PAR     ', myThid)) THEN
        DO l=2,nlam
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx
            PAR(i,j,k,1) = PAR(i,j,k,1) + PAR(i,j,k,l)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        WRITE(diagname, '(A)') 'PAR'
        CALL DIAGNOSTICS_FILL(PAR,diagname,0,Nr,2,bi,bj,myThid)
       ENDIF
#ifndef DARWIN_ALLOW_CHLQUOTA
       IF (DIAGNOSTICS_IS_ON('Chl     ', myThid)) THEN
#ifdef ALLOW_RADTRANS
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx
           tmp3d(i,j,k) = 0 _d 0
          ENDDO
         ENDDO
        ENDDO
        DO l=1,nlam
         DO k=1,Nr
          DO j=1,sNy
           DO i=1,sNx
            tmp3d(i,j,k) = tmp3d(i,j,k) + chlPrev(i,j,k,bi,bj,l)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(tmp3d,'Chl     ',0,Nr,3,bi,bj,myThid)
#else
        CALL DIAGNOSTICS_FILL(chlPrev,'Chl     ',0,Nr,1,bi,bj,myThid)
#endif
       ENDIF
#endif
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iPP),   'PP      ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iNfix), 'Nfix    ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iDenit),'Denit   ',
     &       0,Nr,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iDenitN),'DenitN  ',
     &       0,Nr,2,bi,bj,myThid)

       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          cNO3(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNO3)
          cNO2(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNO2)   
          cNH4(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsNH4)
          cPO4(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsPO4)
          cSi(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsSi)
          cFe(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsFe)
          consDIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsDIC)
          consDIC_PIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsDIC_PIC)
          respirDIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iRespirDIC)
          reminDIC_DOC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iReminDIC_DOC)
          reminDIC_POC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iReminDIC_POC)
          disscDIC_PIC(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iDisscDIC_PIC)    
          cALK(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsALK)
          cO2(i,j,k) = hFacC(i,j,k,bi,bj)
     &     * diags(i,j,k,iConsO2)     
         ENDDO
        ENDDO
       ENDDO

C      nutrient consumption diagnostics
       CALL DIAGNOSTICS_FILL(cNO3(1-OLx,1-OLy,1),
     &  'C_NO3   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cNO2(1-OLx,1-OLy,1),
     &  'C_NO2   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cNH4(1-OLx,1-OLy,1),
     &  'C_NH4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cPO4(1-OLx,1-OLy,1),
     &  'C_PO4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cSi(1-OLx,1-OLy,1),
     &  'C_Si    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cFe(1-OLx,1-OLy,1),
     & 'C_Fe    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(consDIC(1-OLx,1-OLy,1),
     & 'cDIC    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(consDIC_PIC(1-OLx,1-OLy,1),
     & 'cDIC_PIC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(respirDIC(1-OLx,1-OLy,1),
     & 'respDIC ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(reminDIC_DOC(1-OLx,1-OLy,1),
     & 'rDIC_DOC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(reminDIC_POC(1-OLx,1-OLy,1),
     & 'rDIC_POC',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(disscDIC_PIC(1-OLx,1-OLy,1),
     & 'dDIC_PIC',0,Nr,1,bi,bj,myThid)     
       CALL DIAGNOSTICS_FILL(cALK(1-OLx,1-OLy,1),
     & 'C_ALK   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(cO2(1-OLx,1-OLy,1),
     & 'C_O2    ',0,Nr,1,bi,bj,myThid)

C      compute 'source' diagnostics: tendencies without consumption;
C      add full tendencies to consumption and store result back into
C      diags(iCons*) for convenience
       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          sALK(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsALK)
     &     + gPtr(i,j,k,iALK))
          sNO3(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNO3)
     &     + gPtr(i,j,k,iNO3))
          sNO2(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNO2)
     &     + gPtr(i,j,k,iNO2))
          sNH4(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsNH4)
     &     + gPtr(i,j,k,iNH4))
          sPO4(i,j,k) = hFacC(i,j,k,bi,bj) * ((diags(i,j,k,iConsPO4))
     &     + gPtr(i,j,k,iPO4))
          sSi(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsSi)
     &     + gPtr(i,j,k,iSiO2))
          sFe(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsFe)
     &     + gPtr(i,j,k,iFeT))
          sO2(i,j,k) = hFacC(i,j,k,bi,bj) * (diags(i,j,k,iConsO2)
     &     + gPtr(i,j,k,iO2))
 
         ENDDO
        ENDDO
       ENDDO

       CALL DIAGNOSTICS_FILL(sALK(1-OLx,1-OLy,1),
     &  'S_ALK   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNO3(1-OLx,1-OLy,1),
     &  'S_NO3   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNO2(1-OLx,1-OLy,1),
     &  'S_NO2   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sNH4(1-OLx,1-OLy,1),
     &  'S_NH4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sPO4(1-OLx,1-OLy,1),
     &  'S_PO4   ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sSi(1-OLx,1-OLy,1),
     &  'S_Si    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(sFe(1-OLx,1-OLy,1),
     &  'S_Fe    ',0,Nr,1,bi,bj,myThid)
      CALL DIAGNOSTICS_FILL(sO2(1-OLx,1-OLy,1),
     &  'S_O2    ',0,Nr,1,bi,bj,myThid)

      DO iTr=1,nDarwin
       DO k=1,Nr
        DO j=1,sNy
         DO i=1,sNx
          gPtrTemp(i,j,k,iTr) = hFacC(i,j,k,bi,bj)*
     &     gPtr(i,j,k,iTr)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

       DO iTr=1,nDarwin
        diagname = 'gDAR'//PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL(gPtrTemp(1-OLx,1-OLy,1,iTr), diagname,
     &   0,Nr,2,bi,bj,myThid)
       ENDDO

#ifdef DARWIN_DIAG_PERTYPE
       DO iTr=1,nplank
        WRITE(diagname, '(A,I4.4)') 'PP', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iPPplank+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
        WRITE(diagname, '(A,I4.4)') 'PC', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iPCplank+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
        WRITE(diagname, '(A,I4.4)') 'GR', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iGRplank+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
        WRITE(diagname, '(A,I4.4)') 'GrGn', iTr
        CALL DIAGNOSTICS_FILL(diags(1-OLx,1-OLy,1,iGrGn+iTr-1),
     &        diagname,0,Nr,2,bi,bj,myThid)
       ENDDO
#endif
       IF (DIAGNOSTICS_IS_ON('sfcSolFe', myThid)) THEN
        DO j=jMin,jMax
         DO i=iMin,iMax
          IF (hFacC(i,j,1,bi,bj) .EQ. 0.) THEN
           surfFe(i,j) = 0 _d 0
          ELSE
           surfFe(i,j) = alpfe * inputFe(i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO
        CALL DIAGNOSTICS_FILL(surfFe,'sfcSolFe',0,1,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(sedFe,'sedFe   ',0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(diffFe,'freeFe  ',0,Nr,1,bi,bj,myThid)
       ENDIF

#ifdef DARWIN_ALLOW_CARBON
       CALL DIAGNOSTICS_FILL(pH,     'pH      ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(pCO2,   'pCO2    ',0,Nr,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(FluxCO2,'fluxCO2 ',0,1,1,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gDIC,   'gDICsurf',0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gDICsurfForc,'gDICEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gALKsurfForc,'gALKEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gO2surfForc, 'gO2Epr  ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNO3surfForc,'gNO3Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNO2surfForc,'gNO2Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gNH4surfForc,'gNH4Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gPO4surfForc,'gPO4Epr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gFeTsurfForc,'gFeTEpr ',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gSiO2surfForc,'gSiO2Epr',0,1,2,
     &  bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gALK,   'gALKsurf',0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gO2,    'gO2surf ',0,1,2,bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(gO2TH,    'gO2Tsurf',0,1,2,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(gN2O,     'gN2Osurf',0,1,2,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(gN2OTH,   'gN2OTHsu',0,1,2,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(gAR,      'gARsurf ',0,1,2,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(FluxO2TH, 'fluxO2TH',0,1,1,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(FluxN2O,  'fluxN2O ',0,1,1,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(FluxN2OTH,'fluN2OTH',0,1,1,bi,bj,myThid)!makai
       CALL DIAGNOSTICS_FILL(FluxAR,   'fluxAR  ',0,1,1,bi,bj,myThid)!makai
#endif
       CALL TIMER_STOP ('DIAGS_FILL [DARWIN_FORCING]',myThid)
      ENDIF
#endif

      midTime = midTime + deltaTclock/nsubtime
      subTime = subTime + deltaTclock/nsubtime
C     isub
      ENDDO

#endif /* ALLOW_DARWIN */

      RETURN
      END
