!     ==================================================================
!       THIS SUBROUTINE INITIALIZES THE ARRAYS FOR THE PARTICLE PHASE
!     ==================================================================
      
      SUBROUTINE INIT_PPHASE(bMASS_MIN, bMASS_MAX, mSORB, PATHi)
        
      USE hTOMAS; IMPLICIT NONE

!     DECLARATIONS:
!     ==================================================================
!     ABSORBING SEED:
      INTEGER mSORB
      
!     NUMBER OF PARTICLES IN EACH BIN [box-1]:
      DOUBLE PRECISION np
      
!     COUNTERS:
      INTEGER i,j,k

!     MINIMUM [kg] AND
!     MAXIMUM PARTICLE MASS [kg] OF ENTIRE SIZE RANGE:
      DOUBLE PRECISION bMASS_MIN
      DOUBLE PRECISION bMASS_MAX

!     BIN MASS RATIO:
      DOUBLE PRECISION RATIO

!     PATH FOR INPUT FOLDER:
      CHARACTER(LEN=50) :: PATHi
      
!     INITIALIZE THE BOUNDS FOR EACH BIN:
!     ==================================================================
!     BIN MASS RATIO:
      RATIO = EXP(LOG(bMASS_MAX/bMASS_MIN)/DBLE(iBINS))

!     BOUNDS:
      DO k = 1,iBINS+1
         Xk(k) = bMASS_MIN*(RATIO**(k-1))
      END DO
      
!     INITIALIZE THE NUMBER AND MASS DISTRIBUTIONS:
!     ==================================================================
      OPEN(UNIT=12, FILE=TRIM(PATHi)//'/'//'in.PSD0', STATUS='OLD')
      
      DO i = 1,iBINS
         READ(12,*) np

         Nk(i) = &
         Nk(i) + np
         
         IF (mSORB.EQ.0) THEN
            Mk(i,xSO4) = &
            Mk(i,xSO4) + np*SQRT(Xk(i)*Xk(i+1))
         ELSE
            Mk(i,xORG_LST) = &
            Mk(i,xORG_LST) + np*SQRT(Xk(i)*Xk(i+1))
         END IF
      END DO

      CLOSE(12)
      
!     SET TRIVIAL VALUES:
!     ==================================================================
      DO k = 1,iBINS
         IF (Nk(k).LT.1D-60) Nk(k) = 1e-9
         DO j = 1,iCOMP
            IF (Mk(k,j).LT.1D-60) Mk(k,j) = 0.0
         END DO
      END DO

      DO k = 1,iBINS
         DO j = 1,iCOMP-iDIAG
            IF (Mk(k,j).EQ.0.0) THEN
                Mk(k,j) = 1.d-10*SUM(Mk(k,:))
            END IF
         END DO
      END DO
      
      RETURN
      END SUBROUTINE INIT_PPHASE
