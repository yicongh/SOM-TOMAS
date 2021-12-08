!     ======================================================================================
!     THIS SUBROUTINE STEPS FORWARD FOR COAGULATION
!     ======================================================================================
      
      SUBROUTINE STEP_COAG(Temp,Pres,BoxVol,COAG,mdt)
      
      USE hTOMAS; IMPLICIT NONE

!     DECLARATIONS:
!     ======================================================================================
      DOUBLE PRECISION Temp
      DOUBLE PRECISION Pres
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION mdt

!     COAGULATION SWITCH:
      INTEGER COAG

!     PARTICLE DIAMETERS AND
!     DIFFUSIVITIES:
      DOUBLE PRECISION k_DIAM(iBINS)
      DOUBLE PRECISION k_DIFF(iBINS)
            
!     PARTICLE MEAN VELOCITY:
      DOUBLE PRECISION k_MEANV(iBINS)

!     COAGULATION KERNEL:
      DOUBLE PRECISION kk_COAG_KERNEL(iBINS,iBINS)

!     AIR VISCOSITY AND
!     MEAN FREE PATH:
      DOUBLE PRECISION AIR_MU
      DOUBLE PRECISION AIR_MFP

!     KNUDSEN NUMBER:
      DOUBLE PRECISION Kn

!     NON-CONTINUUM CORRECTION FACTOR:
      DOUBLE PRECISION BETA

!     NUMBER OF COAGULATION COLLISIONS:
      DOUBLE PRECISION nnCOAG

!     LOOP COUNTERS:
      INTEGER i,j,k
      
!     INDEX:
      INTEGER INDEX

!     TEMPORARY VARIABLES:
      DOUBLE PRECISION mp1
      DOUBLE PRECISION mp2
      DOUBLE PRECISION mp3

!     BOLTZMANN CONSTANT:
      DOUBLE PRECISION kB
      PARAMETER(kB = 1.38e-23)

!     GAS CONSTANT:
      DOUBLE PRECISION R
      PARAMETER(R = 8.314)

!     PI:
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.141592654)
      
!     SKIPPING CONTROL:
!     ======================================================================================
      IF (COAG.EQ.0) GOTO 101
      
!     CALCULATE THE COAGULATION KERNEL:
!     ======================================================================================
!     AIR VISCOSITY AND
!     MEAN FREE PATH:
      AIR_MU  = 2.5277e-7*(Temp**0.75302)
      AIR_MFP = 2.0*AIR_MU/(Pres*SQRT(8.0*0.0289/(PI*R*Temp)))

!     PARTICLE SIZES AND
!     DIFFUSIVITIES:
      DO k = 1,iBINS
         
         k_DIAM(k) = ((SQRT(Xk(k)*Xk(k+1))/1770.)*(6./PI))**(1./3.)
         
!        KNUDSEN NUMBER:
         Kn = 2.0*AIR_MFP/k_DIAM(k)
         
         k_DIFF(k) = kB*Temp/(3.0*PI*AIR_MU*k_DIAM(k))* &
                    (5.0 + 4.0*Kn + 6.0*Kn**2 + 18.0*Kn**3)/ &
                    (5.0 - Kn + (8.0 + PI)*Kn**2)

!        MEAN VELOCITY:
         k_MEANV(k) = SQRT(8.0*kB*Temp/PI/SQRT(Xk(k)*Xk(k+1)))
         
      END DO

!     COAGULATION KERNEL:
      DO i = 1,iBINS
      DO j = 1,iBINS

!        KNUDSEN NUMBER:
         Kn = 4.0*(k_DIFF(i) + k_DIFF(j))/ &
            SQRT(k_MEANV(i)**2 + k_MEANV(j)**2)/(k_DIAM(i) + k_DIAM(j))

!        NON-CONTINUUM CORRECTION FACTOR:
         BETA = (1.0 + Kn)/(1.0 +2.0*Kn*(1.0 + Kn))
         
!        COAGULATION KERNEL:
         kk_COAG_KERNEL(i,j) = 2.0*PI*(k_DIAM(i) + k_DIAM(j))* &
                               BETA*(k_DIFF(i) +k_DIFF(j))*1e6
         
      END DO
      END DO

!     STEP FOR COAGULATION:
!     ======================================================================================
      DO i = 1,iBINS
      DO j = i,iBINS
         
!        SKIPPING CRITERIA:
         IF (Nk(i).LT.10.) GOTO 200
         IF (Nk(j).LT.10.) GOTO 200
         
!        NUMBER OF COAGULATION EVENTS:
         nnCOAG = kk_COAG_KERNEL(i,j)*Nk(i)/BoxVol*Nk(j)/BoxVol*mdt

!        FIND COMBINED PARTICLE INDEX:
         mp1 = SUM(Mk(i,1:iCOMP-iDIAG))/Nk(i)
         mp2 = SUM(Mk(j,1:iCOMP-iDIAG))/Nk(j)
         mp3 = mp1 + mp2
         
         INDEX = 0
         DO k = j,iBINS
            IF (mp3.LT.Xk(k+1)) THEN
               INDEX = k; GOTO 201
            END IF
         END DO
         INDEX = iBINS
 201     CONTINUE

!        SUBTRACT AND ADD MASS:
         DO k = 1,iCOMP-iDIAG
            
            mp1 = Mk(i,k)/Nk(i)
            mp2 = Mk(j,k)/Nk(j)
            mp3 = mp1 + mp2

            Mk(i,k) = &
            Mk(i,k) - mp1*nnCOAG*BoxVol

            Mk(j,k) = &
            Mk(j,k) - mp2*nnCOAG*BoxVol

            Mk(INDEX,k) = &
            Mk(INDEX,k) + mp3*nnCOAG*BoxVol
            
         END DO
         
!        SUBTRACT AND ADD NUMBER:
         Nk(i) = Nk(i) - nnCOAG*BoxVol
         Nk(j) = Nk(j) - nnCOAG*BoxVol
         
         Nk(INDEX) = & 
         Nk(INDEX) + nnCOAG*BoxVol
         
 200     CONTINUE
      END DO
      END DO

      CALL MNFIX(Nk,Mk)

 101  CONTINUE
      RETURN
      END SUBROUTINE
