C     ==========================================================================
C     THIS FUNCTION STEPS FORWARD FOR AEROSOL-PHASE CHEMISTRY (OLIGOMERIZATION).
C     ==========================================================================
      
      SUBROUTINE STEP_AEROCHEM(k_f,k_r,mOLIG,dt,timer)
      
      USE hTOMAS; IMPLICIT NONE
      
C     DECLERATIONS:
C     ==========================================================================
C     OLIGOMER FORMATION [cm3 molecules-1 s-1] AND
C     DISSOCIATION [s-1] RATES:
      DOUBLE PRECISION k_f
      DOUBLE PRECISION k_r

C     SWITCH FOR OLIGOMERIZATION:
      INTEGER mOLIG

C     TIMESTEP [s] AND
C     TIMER [s]:
      DOUBLE PRECISION dt
      DOUBLE PRECISION timer

C     VOLUME OF CONDENSED PHASE [cm3]:
      DOUBLE PRECISION VolCond
      
C     MOLES OF THE MONOMERS [mol]:
      DOUBLE PRECISION MONO1
      DOUBLE PRECISION MONO2

C     CHANGE IN MOLES OF THE MONOMERS [mol]:
      DOUBLE PRECISION dMOLE1
      DOUBLE PRECISION dMOLE2

C     CHANGE IN MASS OF THE MONOMERS [kg]:
      DOUBLE PRECISION dM1
      DOUBLE PRECISION dM2

C     MOLES OF THE DIMER [mol]:
      DOUBLE PRECISION DIMER

C     CHANGE IN MOLES OF THE DIMER [mol]:
      DOUBLE PRECISION dMOLE

C     CHANGE IN MASS OF THE DIMER [kg]:
      DOUBLE PRECISION dM

C     MOLECULAR WEIGHT OF THE DIMER [g mol-1]:
      DOUBLE PRECISION MW_DIMER

C     AVOGADRO NUMBER [mol-1]:
      DOUBLE PRECISION NA
      PARAMETER(NA=6.022e23)

C     AEROSOL DENSITY [kg m-3]:
      DOUBLE PRECISION RhoORG
      PARAMETER(RhoORG=1400.)
      
C     LOOP COUNTERS:
      INTEGER i,j,k

C     INITIAL MASS:
      DOUBLE PRECISION Mk0(iBINS,iCOMP)
      
C     SKIPPING CRITERIA:
C     ==========================================================================      
      IF (mOLIG.EQ.0) GOTO 100

C     KEEP INITIAL VALUES:
C     ==========================================================================
      DO k = 1,iBINS
         DO j = 1,iCOMP
            Mk0(k,j) = Mk(k,j)
         END DO
      END DO
      
C     STEP FOR THE FORWARD REACTIONS:
C     ==========================================================================
      DO i = 1,iMONO-1
         DO j = i,iMONO-1
            DO k = 1,iBINS
               IF (Nk(k).LT.1.) GOTO 101

C              CALCULATE VOLUME OF CONDENSED PHASE IN BIN:
               VolCond = SUM(Mk(k,xORG_1ST:xORG_LST))/RhoORG*1e6

               IF (VolCond.LT.SQRT(Xk(k)*Xk(k+1))/RhoORG*1e6) GOTO 101
                               
               MONO1 = Mk(k,xORG_1ST+i-1)/(ORG_MW(i)/1000.0)*NA/VolCond
               MONO2 = Mk(k,xORG_1ST+j-1)/(ORG_MW(j)/1000.0)*NA/VolCond
               
               IF (MONO1.LE.0.0) GOTO 101
               IF (MONO2.LE.0.0) GOTO 101

               dMOLE1 = MONO1*(1.0 - EXP(-k_f*MONO2*dt))
               dMOLE2 = MONO2*(1.0 - EXP(-k_f*MONO1*dt))
               
               dM1 = dMOLE1*VolCond/NA*(ORG_MW(i)/1000.0)
               dM2 = dMOLE2*VolCond/NA*(ORG_MW(j)/1000.0)
                              
C              ADD TO DIMERS:
               Mk(k,xDIME_1ST+i-1) = Mk(k,xDIME_1ST+i-1) + dM1
               Mk(k,xDIME_1ST+j-1) = Mk(k,xDIME_1ST+j-1) + dM2

C              DEDUCT FROM MONOMERS:
               Mk(k,xORG_1ST+i-1) = Mk(k,xORG_1ST+i-1) - dM1
               Mk(k,xORG_1ST+j-1) = Mk(k,xORG_1ST+j-1) - dM2
 101           CONTINUE
            END DO
         END DO
      END DO

C     STEP FOR THE REVERSE REACTIONS:
C     ==========================================================================
      DO i = 1,iMONO-1
         DO k = 1,iBINS
            DIMER = Mk(k,xDIME_1ST+i-1)/(ORG_MW(i)/1000.)
            dMOLE = DIMER*(1.0 - EXP(-k_r*dt))
            dM    = dMOLE*(ORG_MW(i)/1000.0)

            Mk(k,xDIME_1ST+i-1) = Mk(k,xDIME_1ST+i-1) - dM
            Mk(k,xORG_1ST+i-1)  = Mk(k,xORG_1ST+i-1)  + dM            
         END DO
      END DO

C     CALCULATE OLIG. CONTRIBUTION TO PARTICLE-PHASE DECAY:
C     ==========================================================================
C     REFRESH THE ARRAY:
      kc_OLIG(:,:) = 0.0

C     PSEUDO-FIRST-ORDER DECAY:         
      DO k = 1,iBINS
         DO j = 1,iCOMP
            kc_OLIG(k,j) = LOG(MIN(Mk(k,j)/Mk0(k,j),1.0))/(-dt)
         END DO
      END DO

C     FIX SMALL VALUES:
C     ==========================================================================
      DO k = 1,iBINS
         IF (Nk(k).LE.1D-5) THEN
            Nk(k) = 1D-5
         END IF
         DO j = 1,iCOMP
            IF (Mk(k,j).LE.1D-5*SQRT(Xk(k)*Xk(k+1))) THEN 
               Mk(k,j) = 1D-5*SQRT(Xk(k)*Xk(k+1))
            END IF
         END DO
      END DO
      
 100  RETURN
      END SUBROUTINE
