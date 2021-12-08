C     ===================================================================
C     THIS FUNCTION CALCULATES THE AEROSOL-SIDE MASS TRANSFER COEFFICIENT
C     ===================================================================
      
      DOUBLE PRECISION FUNCTION CALC_KPAR(Rc,Rp,Db,kc)
      
      USE hBTAB
      
      DOUBLE PRECISION Rc
      DOUBLE PRECISION Rp
      DOUBLE PRECISION Db
      DOUBLE PRECISION kc
      DOUBLE PRECISION QFACTOR
      
      DOUBLE PRECISION OUT1
      DOUBLE PRECISION QBULK
      
      DOUBLE PRECISION OUT2
      DOUBLE PRECISION KPAR

      DOUBLE PRECISION L
      DOUBLE PRECISION BETAS(10)
      
      INTEGER i,j,k
      INTEGER INDEX1,INDEX2
      
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.1415926)
      
      DOUBLE PRECISION X1,X2
      DOUBLE PRECISION Y1(10)
      DOUBLE PRECISION Y2(10)
      
      DOUBLE PRECISION SIGMA
      DOUBLE PRECISION G

C     Q-FACTOR:
      QFACTOR = (Rp - Rc)*SQRT(kc/Db)
      
C     CALCULATE THE L VALUE:
      L = (Rp - Rc)/Rc

C     IF L IS GREATER THAN 1000:
      IF (L.GE.1000.0) THEN
         DO i = 1,10
            BETAS(i) = DBLE(i)*PI
         END DO
         GOTO 100
      END IF

C     OTHERWISE INTERPOLATE TO FIND BETAS:
      INDEX1 = 1

      DO i = 1,61
         IF (L_TAB(INDEX1).LE.L) THEN
            INDEX1 = INDEX1 + 1
         ELSE
            INDEX1 = i - 1
            GOTO 200
         END IF
      END DO
      
 200  CONTINUE
      INDEX2 = INDEX1 + 1

      X1 = L_TAB(INDEX1)
      X2 = L_TAB(INDEX2)

      Y1 = BETAS_TAB(INDEX1,:)
      Y2 = BETAS_TAB(INDEX2,:)

      BETAS = Y1 + (Y2 - Y1)/(X2 - X1)*(L - X1)
      
 100  CONTINUE
      
C     CALCULATE THE Q-BULK FACTOR:
      SIGMA = 0.0
      
      DO i = 1,10
         G = Rc*BETAS(i)*COS(BETAS(i)) + (Rp-Rc)*SIN(BETAS(i)) -
     &       Rp*BETAS(i)
         
         SIGMA = SIGMA + (BETAS(i)**2.0 + L**2.0)*G/
     &                   (L + L**2.0 + BETAS(i)**2.0)/(BETAS(i)**3.0)/
     &                   (BETAS(i)**2.0 + QFACTOR**2.0)
      END DO

      QBULK = 1.0 + 6.0*Rp*(Rp - Rc)/
     &       (Rp**3.0 - Rc**3.0)*SIGMA*(QFACTOR**2.0)

      KPAR = Db*QBULK/
     &      (-6.0*Rp*((Rp - Rc)**3.0)*SIGMA/(Rp**3.0 - Rc**3.0))*
     &      (Rp**3.0 - Rc**3.0)/3.0/Rp**2.0
      
      CALC_KPAR = KPAR

      RETURN
      END FUNCTION
