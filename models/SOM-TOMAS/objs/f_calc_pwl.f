C ====================================================================
C This function calculates the size-dependent particle wall loss rates
C ====================================================================

      DOUBLE PRECISION FUNCTION CALC_PWL(Dp)

C     DECLARATIONS:
C     ================================================================
      IMPLICIT NONE
      
C     SIZE OF PARTICLE [nm],
C     CUTOFF SIZE FOR FIT [nm]:
      DOUBLE PRECISION Dp
      DOUBLE PRECISION CUT
      PARAMETER(CUT = LOG10(210.))

C     POLYNOMIAL PARAMETERIZATIONS:
      DOUBLE PRECISION a1,a2,a3,a4
      PARAMETER(a1 =  0.02727767)
      PARAMETER(a2 =  0.25517087)
      PARAMETER(a3 = -3.04476291)
      PARAMETER(a4 =  1.44203902)
      
      DOUBLE PRECISION b1,b2,b3,b4
      PARAMETER(b1 =  -13.3175139)
      PARAMETER(b2 =  103.60122852)
      PARAMETER(b3 = -264.61991781)
      PARAMETER(b4 =  218.59323301)
      
C     WALL LOSS RATE [s-1]:
      DOUBLE PRECISION k_par

C     FUNCTION:
C     ================================================================
C     CONVERT SIZE TO LOG SPACE:
      IF (Dp.LT.20.0) Dp = 20.0
      
      Dp = LOG10(Dp)
      
C     FIND LOSS RATE:
      IF (Dp.LE.CUT) THEN
         k_par = a1*Dp**3. + a2*Dp**2. + a3*Dp  + a4
         GOTO 10
      END IF
      
      IF (Dp.GT.CUT) THEN
         k_par = b1*Dp**3. + b2*Dp**2. + b3*Dp  + b4
         GOTO 10
      END IF
      
 10   CONTINUE

C     CONVERT LOSS RATE BACK TO LINEAR SPACE:
      CALC_PWL = (10.0**(k_par))/60.
      
      RETURN
      END FUNCTION
