!     ============================================================================
!                       THIS SUBROUTINE WRITES THE OUTPUT FILES
!     ============================================================================

      SUBROUTINE WRITE_OUTPUTS(timer,Db,OH,BoxVol,JNUC)

      USE hSOM; USE hTOMAS; USE hOFILES; IMPLICIT NONE
     
!     DELARATIONS:
!     ============================================================================
      DOUBLE PRECISION timer
      DOUBLE PRECISION Db
      DOUBLE PRECISION OH
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION JNUC
      
!     COUNTERS:
      INTEGER i,j,k

!     VARIABLE NAMES:
      DOUBLE PRECISION O2C,OM2OC,CNO,ONO
      DOUBLE PRECISION SOA1,SOA2,SOA3,OLIG

      DOUBLE PRECISION NUM1,NUM2
      DOUBLE PRECISION dVOC

!     FORMATS:
      CHARACTER*120 FMT1,FMT2,FMT3,FMT4,FMT5,FMT6,FMT7
      
      DOUBLE PRECISION ELVOC

!     SPECIFY THE FORMATS:
!     ============================================================================
      FMT1 = '(0007E15.6)'
      FMT2 = '(0003E15.6)'
      FMT3 = '(0002E15.6)'
      FMT4 = '(0000E15.6)'; WRITE(FMT4(2:5),'(I4)') iBINS + 1 
      FMT5 = '(0000E15.6)'; WRITE(FMT5(2:5),'(I4)') nSOM + 1
      FMT6 = '(0000E15.6)'; WRITE(FMT6(2:5),'(I4)') iCOMP + 1
      FMT7 = '(0000E15.6)'; WRITE(FMT7(2:5),'(I4)') iCOMP - 1 + 1
      
!     SET TRIVIAL VALUES:
!     ============================================================================
      DO i = 1,iBINS
         IF (Nk(i).LT.1d-60) Nk(i) = 1d-60
         DO j = 1,iCOMP
            IF (Mk(i,j).LT.1d-60) Mk(i,j) = 1d-60
         END DO
      END DO

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_SOA,STATUS='OLD',ACCESS='APPEND')
      
      SOA1 = SUM(Mk(:,xORG_1ST:xORG_LST-1))*1e9/(BoxVol*1e-6)      
      SOA2 = SUM(Mk_WALL(:,xORG_1ST:xORG_LST))*1e9/(BoxVol*1e-6)

      OLIG = SUM(Mk(:,xDIME_1ST:xORG_LST-1))*1e9/(BoxVol*1e-6)
      
      ONO = SUM(SUM(Mk(:,xORG_1ST:xORG_LST-1),1)/ORG_MW(1:iORG-1)*ONUMBER(1:iORG-1))
      CNO = SUM(SUM(Mk(:,xORG_1ST:xORG_LST-1),1)/ORG_MW(1:iORG-1)*CNUMBER(1:iORG-1))
      
      O2C = ONO/CNO
      
      WRITE(35,FMT1) timer/3600.,SOA1,SOA2,O2C,ONO,CNO,OLIG
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_NCONC,STATUS='OLD',ACCESS='APPEND')
      
      NUM1 = SUM(Nk(:))/BoxVol
      NUM2 = SUM(Nk_wall(:))/BoxVol
      
      WRITE(35,FMT2) timer/3600.,NUM1,NUM2
      CLOSE(35)
      
!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_PREC,STATUS='OLD',ACCESS='APPEND')
      WRITE(35,FMT3) timer/3600.,SOMGC(xPREC_1ST(1))
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_MDIST,STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(35,FMT4) timer/3600., &
           (SUM(Mk(k,1:xORG_LST))*1e9/BoxVol/1e-6/ &
            LOG10(Xk(k+1)/Xk(k))/3.0,k=1,iBINS)
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_NDIST,STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(35,FMT4) timer/3600., &
           (Nk(k)/BoxVol/(LOG10(Xk(k+1)/Xk(k))/3.0),k=1,iBINS)
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_SOMGC,STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(35,FMT5) timer/3600., (SOMGC(j),j=1,nSOM)
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_DB,STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(35,FMT3) timer/3600., DB
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE='outputs/out.JNUC',STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(35,FMT3) timer/3600., MAX(JNUC,1D-60)
      CLOSE(35)
      
!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_AEMASS,STATUS='OLD',ACCESS='APPEND')
      WRITE(35,FMT6) timer/3600.,(SUM(Mk(:,j))*1e9*1e6/BoxVol,j=1,iCOMP)
      CLOSE(35)

!     SAVE DATA:
!     ============================================================================
      OPEN(UNIT=35,FILE=OUT_GCMASS,STATUS='OLD',ACCESS='APPEND')
      WRITE(35,FMT7) timer/3600.,(Gc(j)*1e9*1e6/BoxVol,j=1,iCOMP-1)
      CLOSE(35)
      
      RETURN
      END SUBROUTINE
