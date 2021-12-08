C     =========================================================================
C         THIS SUBROUTINE SOLVES SIZE DISTRIBUTION NUMBER AND MASS MOMENTS
C     =========================================================================
      
      SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,xSPECIES,MASK)
      
      USE hTOMAS; IMPLICIT NONE
      
C     DECLARATIONS:
C     =========================================================================
C     GROWTH FORCE AND BIN BOUNDS:
      DOUBLE PRECISION TAU(iBINS)
      DOUBLE PRECISION X(iBINS+1)

C     NUMBER AND MASS ARRAYS:
      DOUBLE PRECISION ANK(iBINS)
      DOUBLE PRECISION ANKD(iBINS)
      DOUBLE PRECISION AMK(iBINS,iCOMP)
      DOUBLE PRECISION AMKD(iBINS,iCOMP)
      
      DOUBLE PRECISION AMKDRY(iBINS)
      DOUBLE PRECISION AMKWET(iBINS)
      DOUBLE PRECISION WR(iBINS)

C     INDEX OF CONDENSING SPECIES:
      INTEGER xSPECIES

C     CONDENSATION MASK:
      INTEGER MASK(iBINS)

C     COUNTERS:
      INTEGER i,j,k
      INTEGER IMN

C     TRAZOID APPROX. VARIABLES:
      DOUBLE PRECISION EFF,PSI
      DOUBLE PRECISION RATIO
      DOUBLE PRECISION AVG
      DOUBLE PRECISION DIV

C     BIN MASS VARIABLES:
      DOUBLE PRECISION DN,DM
      DOUBLE PRECISION XL,XU
      DOUBLE PRECISION YL,YU
      
      INTEGER INDEX1,INDEX2
   
      DOUBLE PRECISION AA,A1,A2,HH
      DOUBLE PRECISION MM,M1,M2
      DOUBLE PRECISION dMASS
      
C     MAX OF TAU:
      DOUBLE PRECISION MAXTAU
      DOUBLE PRECISION,EXTERNAL :: GET_MAX
      
C     MINIMUM TAU:
      DOUBLE PRECISION Teps
      PARAMETER(Teps = 1.0d-40)

C     MINIMUM NUMBER IN BIN:
      DOUBLE PRECISION Neps
      PARAMETER(Neps = 1.0d-60)

C     MASS GROWTH FUNCTION:
      DOUBLE PRECISION,EXTERNAL :: DMDT_INT

      DOUBLE PRECISION VV1,VV2      
      DOUBLE PRECISION XX1,XX2
      
C     INITIALIZE ARRAYS:
C     =========================================================================
      DO k = 1,iBINS
         IF (ANKD(k).LT.Neps) THEN
            ANKD(k) = Neps
            AMKD(k,xSO4) = Neps*SQRT(X(k)*X(k+1))
            DO j = 1,iCOMP
               IF (j.NE.xSO4) THEN
                  AMKD(k,j) = 1D-60
               END IF
            END DO
         END IF
      END DO
      
      DO k = 1,iBINS
         AMKDRY(k)=0.0d0
         AMKWET(k)=0.0d0
         DO j=1,iCOMP-iDIAG
            AMKDRY(k)=AMKDRY(k)+AMKD(k,j)
         ENDDO
         DO j=1,iCOMP
            AMKWET(k)=AMKWET(k)+AMKD(k,j)
         ENDDO
         WR(k)=AMKWET(k)/AMKDRY(k)
      ENDDO 
      
      DO k = 1,iBINS
         DO j = 1,iCOMP
            AMK(k,j) = 1.d-60
         END DO
         ANK(k) = 1.d-60
      END DO
      
C     LOOP OVER BINS:
C     =========================================================================
      MAXTAU = GET_MAX(iBINS,TAU)
      
      IF (ABS(MAXTAU).LT.TEPS) THEN
         DO i = 1,iBINS
            DO j = 1,iCOMP
               AMK(i,j) = AMKD(i,j)
            END DO
            ANK(i) = ANKD(i)
         END DO   
         GOTO 300
      END IF
      
      DO 200 k = 1,iBINS
      
      IF (MASK(k).EQ.0) THEN
         ANK(k) = ANK(k) + ANKD(k)
         DO j = 1,iCOMP
            AMK(k,j) = AMK(k,j) + AMKD(k,j)
         ENDDO
         GOTO 200
      END IF

      IF (TAU(k).EQ.0.) THEN
         ANK(k) = ANK(k) + ANKD(k)
         DO j = 1,iCOMP
            AMK(k,j) = AMK(k,j) + AMKD(k,j)
         END DO
         GOTO 200
      END IF

C     TRAPEZOID APPROX. COEFFICIENTS:
C     =========================================================================
      RATIO = X(2)/X(1)
      
      AVG = AMKDRY(k)/ANKD(k)
      DIV = AVG/X(k)
         
      IF (DIV.GT.RATIO) DIV = RATIO
      IF (DIV.LT.1.0  ) DIV = 1.0
         
      EFF = 2.0*ANKD(k)/X(k)*(RATIO-DIV)*(1.0/(RATIO-1.0)**2.0)
      PSI = 2.0*ANKD(k)/X(k)*(DIV  -1.0)*(1.0/(RATIO-1.0)**2.0)

      IF (ANKD(k).LT.1.0) THEN
         ANK(k) = ANK(k) + ANKD(k)
         DO j = 1,iCOMP
            AMK(k,j) = AMK(k,j) + AMKD(k,j)
         END DO
         GOTO 200
      END IF
         
C     GROW THE BIN BOUNDS:
C     =========================================================================
      XU = X(k+1)
      XL = X(k)

      YU = DMDT_INT(XU,TAU(k),WR(k))
      YL = DMDT_INT(XL,TAU(k),WR(k))

      IF (YU.GT.X(iBINS+1)) THEN
         YU  = X(iBINS+1)
      END IF
      
C     BEGIN REMAPPING:
C     =========================================================================
      dMASS = (EFF*YL + PSI*YU)*(YU - YL)*0.5 -
     &        (EFF*XL + PSI*XU)*(XU - XL)*0.5
      
      IF (YU.LT.X(1)) THEN
         ANK(k) = ANK(k) + ANKD(k)
         DO j = 1,iCOMP
            AMK(k,j) = AMK(k,j) + AMKD(k,j)
         END DO
         GOTO 200
      END IF
         
      IF (YL.LT.X(1)) THEN
         ANK(k) = ANK(k) + ANKD(k)
         DO j = 1,iCOMP
            IF (j.EQ.xSPECIES) THEN
               AMK(k,j) = AMKD(k,j) + dMASS
            ELSE
               AMK(k,j) = AMK(k,j) + AMKD(k,j)
            END IF
         END DO
         GOTO 200   
      END IF
         
C     BIN INDEXES:
      IF (TAU(k).GT.0.) THEN
         IMN = k
      ELSE
         IMN = 1
      END IF
      
      DO I = IMN,iBINS
         IF (YL.LE.X(I+1)) THEN
            INDEX1 = I
            GOTO 202
         END IF
      END DO
 202  CONTINUE
      INDEX2 = INDEX1 + 1

C     TRAPEZOIDAL AREAS:
      AA = (EFF + PSI)*(YU - YL)*0.5
      HH = EFF + (PSI - EFF)*(X(INDEX1+1) - YL)/(YU - YL)
      A1 = 0.5*(HH + EFF)*(X(INDEX1+1)-YL)
      A2 = AA - A1
      
      MM = (EFF*YL + PSI*YU)*(YU - YL)*0.5
      HH = (EFF*YL + (PSI*YU - EFF*YL)*(X(INDEX1+1) - YL)/(YU - YL))
      M1 = 0.5*(HH + EFF*YL)*(X(INDEX1+1)-YL)
      M2 = MM - M1

      XX1 = M1/MM
      XX2 = M2/MM
      
      IF (XX1.LT.0.0) XX1 = 0.0 
      IF (XX2.LT.0.0) XX2 = 0.0

      IF (XX1.GT.1.0) XX1 = 1.0      
      IF (XX2.GT.1.0) XX2 = 1.0

      VV1 = A1/AA
      VV2 = A2/AA
      
      IF (VV1.LT.0.0) VV1 = 0.0 
      IF (VV2.LT.0.0) VV2 = 0.0

      IF (VV1.GT.1.0) VV1 = 1.0      
      IF (VV2.GT.1.0) VV2 = 1.0
     
C     CURRENT BIN:
C     =========================================================================
      ANK(INDEX1) = ANK(INDEX1) + ANKD(k)*VV1
      
      DO j = 1,iCOMP
         IF (j.EQ.xSPECIES) THEN
            AMK(INDEX1,j) = 
     &      AMK(INDEX1,j) + (AMKD(k,j) + dMASS)*XX1
         ELSE
            AMK(INDEX1,j) = 
     &      AMK(INDEX1,j) + AMKD(k,j)*XX1
         END IF
      END DO
      
C     NEXT BIN:
C     =========================================================================
      ANK(INDEX2) = ANK(INDEX2) + ANKD(k)*VV2
      
      DO j = 1,iCOMP
         IF (j.EQ.xSPECIES) THEN
            AMK(INDEX2,j) = 
     &      AMK(INDEX2,j) + (AMKD(k,j) + dMASS)*XX2
         ELSE
            AMK(INDEX2,j) = 
     &      AMK(INDEX2,j) + AMKD(k,j)*XX2
         ENDIF
      END DO
      
200   CONTINUE
            
      DO k = 1,iBINS
         IF (ANK(k).LE.1D-10) THEN
            ANK(k) = 1D-10
         END IF
         DO j = 1,iCOMP
            IF (AMK(k,j).LE.1D-10*SQRT(X(k)*X(k+1))) THEN 
               AMK(k,j) = 1D-10*SQRT(X(k)*X(k+1))
            END IF
         END DO
      END DO
      
300   RETURN
      END
