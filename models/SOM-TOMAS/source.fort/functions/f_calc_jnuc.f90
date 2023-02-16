!     ========================================================================================
!                           THIS FUNCTION CALCULATES THE NUCLEATION RATE
!     ========================================================================================

      DOUBLE PRECISION FUNCTION CALC_JNUC(NUCL_A,NUCL_B,NUCL_C,Pres,Temp)
      
      USE hSOM; USE hTOMAS; IMPLICIT NONE      
      
!     NUCLEATION PARAMETERIZATIONS:
      DOUBLE PRECISION NUCL_A
      DOUBLE PRECISION NUCL_B
      DOUBLE PRECISION NUCL_C

!     PRESSURE AND TEMPERATURE:
      DOUBLE PRECISION Pres
      DOUBLE PRECISION Temp

!     MAIN-PROGRAM ITERATION COUNTER:
      INTEGER iCTR

!     LOOP COUNTERS:
      INTEGER i,j,k

!     INDICES:
      INTEGER INDEX,INDEX1,INDEX2
      
!     ELVOC CONCENTRATION:
      DOUBLE PRECISION ELVOC
      
!     CALCULATE ELVOC CONCENTRATION:
!     ========================================================================================
      ELVOC = 0.0
      
      DO i = 1,iMONO
!         IF (CSAT_REF(i).LT.1e-2) THEN
!         IF ((i.GE.iMONO-20).AND.(i.LE.iMONO-14)) THEN
         IF (i.EQ.iMONO-1) THEN            
            
            INDEX = xGRID_1ST(1) + i - 1
            
!            PRINT*, SOMNAMES(INDEX)
            
            ELVOC = &
            ELVOC + SOMGC(INDEX)*1e-6*Pres/Temp/8.314*1e-6*6.022e23
            
         END IF
      END DO
      
!      STOP
      
!     CALCULATE NUCLEATION RATE:
!     ========================================================================================
      IF (ELVOC.LT.10) THEN
         CALC_JNUC = 0.0
      ELSE
         CALC_JNUC = NUCL_A*((ELVOC/NUCL_C)**NUCL_B)
      END IF
      
      RETURN
      END FUNCTION CALC_JNUC
