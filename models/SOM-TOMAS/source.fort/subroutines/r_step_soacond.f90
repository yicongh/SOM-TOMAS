!     ====================================================================================
!     This subroutine determines the condensational driving force of organics
!     between gas and aerosol phases. It then calls the algorithm for
!     condensation (SOLVE_COND).
!      
!     WRITTEN BY Jeff Pierce, Sep. 2010
!     ====================================================================================
      
      SUBROUTINE STEP_SOACOND(Pres,Temp,RH,BoxVol,Alpha,Db,mdt)

      USE hTOMAS; USE hSOM; IMPLICIT NONE

      DOUBLE PRECISION Pres
      DOUBLE PRECISION Temp
      DOUBLE PRECISION RH
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION Alpha
      DOUBLE PRECISION Db
      
      DOUBLE PRECISION mdt
      DOUBLE PRECISION mmdt
      INTEGER,PARAMETER :: nSTEP = 1
      
      DOUBLE PRECISION Nki(iBINS), Mki(iBINS,iCOMP), Gci(iCOMP-1)
      DOUBLE PRECISION Nkf(iBINS), Mkf(iBINS,iCOMP), Gcf(iCOMP-1)
      DOUBLE PRECISION Nko(iBINS), Mko(iBINS,iCOMP), Gco(iCOMP-iDIAG)
      
      DOUBLE PRECISION,PARAMETER :: ST_ORG = 0.025

      DOUBLE PRECISION MS      
      DOUBLE PRECISION MFP
      DOUBLE PRECISION FUCHS
      DOUBLE PRECISION KELVIN
      DOUBLE PRECISION Kn
      
      DOUBLE PRECISION kc      
      DOUBLE PRECISION DGAS
      
      DOUBLE PRECISION KGAS
      DOUBLE PRECISION KTOT
      DOUBLE PRECISION KPAR
      DOUBLE PRECISION QBULK(iBINS)

      DOUBLE PRECISION,EXTERNAL :: CALC_KPAR
      DOUBLE PRECISION,EXTERNAL :: CALC_QBULK
      
      DOUBLE PRECISION pAMB
      DOUBLE PRECISION pSAT
      DOUBLE PRECISION SRATIO      
      DOUBLE PRECISION dPRES(ibins)

      DOUBLE PRECISION ATAU      
      DOUBLE PRECISION TAUj
      DOUBLE PRECISION TAUk(iBINS)
      DOUBLE PRECISION TAU(iBINS)

      DOUBLE PRECISION PSAT_ORG(iMONO)
      DOUBLE PRECISION CSAT_ORG(iMONO)

      DOUBLE PRECISION DENSITY      
      DOUBLE PRECISION DIAM(ibins)
      DOUBLE PRECISION Rp
      DOUBLE PRECISION Rc
      INTEGER,PARAMETER :: CORE_SHELL = 1
      
      INTEGER MASK(iBINS)
      
      DOUBLE PRECISION CONDSINK
      DOUBLE PRECISION SORBPHASE
      DOUBLE PRECISION,PARAMETER :: NONORGSCALE = 0.0
      
      DOUBLE PRECISION mc
      DOUBLE PRECISION dMASS
      DOUBLE PRECISION dMASS2EQ(iBINS)
      
      DOUBLE PRECISION MASS_TOT
      DOUBLE PRECISION MASS_AMS
      DOUBLE PRECISION MASS_ORG
      DOUBLE PRECISION pMASS
      
      INTEGER i,k,j,jj,jo
      INTEGER CTR

      DOUBLE PRECISION,PARAMETER :: PI = 3.141592654
      DOUBLE PRECISION,PARAMETER :: Na = 6.023e23
      DOUBLE PRECISION,PARAMETER :: R = 8.314
      
!     VOLATILITIES:
!     ====================================================================================      
      DO j = 1,iMONO
         PSAT_ORG(j) = PSAT_REF(j)*EXP(Hvap(j)*1.e3*(1.0/298.0 - 1.0/Temp)/R)
         CSAT_ORG(j) = CSAT_REF(j)
      END DO
      
!     SWAP ARRAYS:
!     ====================================================================================
      DO k = 1,iBINS
         Nkf(k) = Nk(k)
         DO j = 1,iCOMP
            Mkf(k,j) = Mk(k,j)
         END DO
      END DO
      
      DO j = 1,iCOMP-1
         Gcf(j) = Gc(j)
      END DO
      
!     LOOP OVER EACH SPECIES:
!     ====================================================================================
      LOOP1: DO jo = 1,iMONO

      j = xORG_1ST + jo - 1
      
!     MEAN SPEED:
      MS = SQRT(8.0*R*Temp/(PI*ORG_MW(jo)*1.0e-3))

!     GAS-PHASE DIFFUSIVITY: 
      DGAS = 1.38e-5*44.0/ORG_MW(jo)

!     MEAN FREE PATH:      
      MFP = 3.d0*DGAS/MS

!     TAU FACTOR:      
      TAUj = 2.0*PI*ORG_MW(jo)*1.0e-3/(R*Temp)
      
!     LOOP OVER INTERNAL TIMESTEPS:
!     ====================================================================================      
      CTR = 0
      
100   CONTINUE

      CALL MNFIX(Nkf,Mkf)
      
!     COND. DRIVING FORCE FACTOR AND TOTAL CONDENSATION SINK:
!     ====================================================================================
      CONDSINK = 0.0

      DO k = 1,iBINS
!        DENSITY OF PARTICLES AND
!        MASS PER PARTICLE:
         IF (Nkf(k).GT.1e-3) THEN
            MASS_ORG = SUM(Mk(k,xORG_1ST:xORG_LST))
            MASS_AMS = Mk(k,xSO4)
            MASS_TOT = MASS_ORG + MASS_AMS
            
            DENSITY = MASS_TOT/(MASS_ORG/1400. + MASS_AMS/1770.); pMASS = MASS_TOT/Nkf(k)
         ELSE
            DENSITY = 1700.; pMASS = SQRT(Xk(k)*Xk(k+1))
         END IF
            
!        PARTICLE DIAMETER:         
         DIAM(k)=((pMASS/DENSITY)*(6./PI))**(0.3333)

!        PARTICLE AND CORE RADII:
         Rp = 0.5*DIAM(k)
         Rc = MIN(0.9999*Rp,0.5*((Mk(k,xSO4)/Nkf(k)/DENSITY)*(6./PI))**(0.3333))

!        KNUDSEN NUMBER:
         Kn = MFP/Rp

!        FUCHS-SUTUGIN CORRECTION:
         FUCHS = 0.75*Alpha*(1.0 + Kn)/(Kn*(1.0 + Kn) + 0.238*Alpha*Kn + 0.75*Alpha)
         
!        FIRST-ORDER DECAY RATE:         
         kc = kc_HET(k,xORG_1ST+jo-1) + kc_OLIG(k,xORG_1ST+jo-1)
         
!        GAS-SIDE MASS-TRANSFER COEFFICIENT:         
         KGAS = DGAS/Rp*FUCHS
         
!        PARTICLE-SIDE MASS-TRANSFER COEFFICIENT: 
         KPAR = CALC_KPAR(Rc,Rp,Db,kc)

!        BULK-TO-SURFACE RATIO:         
         QBULK(k) = CALC_QBULK(Rc,Rp,Db,kc)
                  
!        OVERALL MASS-TRANSFER COEFFICIENT:
         KTOT = 1.0/(1.0/KGAS + 1.0/KPAR*CSAT_ORG(jo)/(DENSITY*1.0e9)*QBULK(k))
         
!        TAU FACTOR:         
         TAUk(k)=(6./(PI*DENSITY))**(1./3.)*KTOT*Rp

!        CUMULATIVELY ADD TO CONDENSATION SINK:         
         IF (Nkf(k).GT.0.0) THEN
            CONDSINK = CONDSINK + PI*(DIAM(k)**2.0)*Nkf(k)*KTOT/(BoxVol*1e-6)
         END IF
         
      END DO

!     PRESSURE DIFFERENCE BETWEEN GAS-PHASE AND PARTICLE SURFACE:
!     ====================================================================================
      DO k = 1,iBINS
!        ABSORBING PHASE:
         SORBPHASE = SUM(Mk(k,xORG_1ST:xORG_LST)) + Mk(k,xSO4)*NONORGSCALE

!        KELVIN RATIO:         
         KELVIN = EXP((ST_ORG*0.001*ORG_MW(jo))/(8.314*Temp*1200.*DIAM(k)))

!        GAS-PHASE PRESSURE:         
         pAMB = (Gcf(xORG_1ST+jo-1)/ORG_MW(jo))/(0.0289*Pres*BoxVol*1.0e-6/R/Temp/28.9)*Pres 
         
!        PARTICLE SURFACE PRESSURE:         
         IF (SORBPHASE.GT.0.0) THEN
            pSAT = PSAT_ORG(jo)*Mkf(k,xORG_1ST+jo-1)/SORBPHASE*KELVIN/QBULK(k)
         ELSE
            pSAT = 0.0
         END IF

!        PRESSURE DIFFERENCE:         
         dPRES(k) = pAMB - pSAT
         
!        MASS CHANGE TO REACH EQUILIBRIUM: 
         IF (SORBPHASE.EQ.0.0) THEN
            dMASS2EQ(k) = 0.0
         ELSE
            dMASS2EQ(k) = (pAMB*QBULK(k)/(PSAT_ORG(jo)/SORBPHASE*KELVIN) - Mkf(k,xORG_1ST+jo-1))/Nkf(k)
         END IF
       
      END DO

!     DETERMINE THE INTEGRATION TIMESTEP:
!     ====================================================================================
      mmdt = mdt/DBLE(nSTEP)

!     CALCULATE COND. DRIVING FORCE:
!     ====================================================================================
      TAU(:) = 0.0
      
      DO k = 1,iBINS
!        CHECK FOR NAN: 
         IF (ISNAN(CONDSINK)) THEN
            PRINT*,'NAN in CONDSINK for j = ', jo; STOP
         END IF
         
         IF (CONDSINK.GT.0.0) THEN
            ATAU = TAUj*TAUk(k)*dPRES(k)/CONDSINK*(1.0 - EXP(-1.0*CONDSINK*mmdt))
         ELSE
            ATAU = 0.0
         END IF
         
!        REGULATE TAU TO EQUILIBRIUM:
         mc = SUM(Mkf(k,:))/Nkf(k)
         
         IF ((ATAU/1.5 + mc**(2.d0/3.d0)).GE.0.0) THEN
            dMASS = (ATAU/1.5 + mc**(2.d0/3.d0))**1.5 - mc
         ELSE
            dMASS = dMASS2EQ(k)
         END IF

         IF ((dMASS/dMASS2EQ(k)).GT.1e-30) THEN
            dMASS = dMASS2EQ(k)*(1.0 - EXP(-dMASS/dMASS2EQ(k)))
         ELSE
            dMASS = 0.0
         END IF
         
         IF (Nkf(k).LT.1e-5) THEN
            TAU(k) = 0.0
         ELSE
            TAU(k) = 1.5*((mc + dMASS)**(2.d0/3.d0) - mc**(2.d0/3.d0))
         END IF
         
      END DO
      
!     INTEGRATION FOR GAS-PARTICLE PARTITIONING:
!     ====================================================================================
!     CONDENSATION MASK:
      MASK(:) = 1

!     INTEGRATION:      
      CALL TMCOND(TAU,Xk,Mkf,Nkf,Mko,Nko,j,MASK)

!     NAN CHECK:
      DO k = 1,iBINS
         IF (ISNAN(Nko(k))) THEN
            PRINT*, 'NaN in Nk after TMCOND'; STOP
         END IF
      END DO

!     CONSERVATION CHECK:
      IF (ABS(SUM(Nko)-SUM(Nkf))/SUM(Nkf).GT.1e-4) THEN
         PRINT*, 'Number Not Conserved in TMCOND'; STOP
      END IF

!     NEGATIVE CHECK:      
      DO k = 1,iBINS
         IF (Nko(k).LT.0.0) THEN
            PRINT*, 'Negative Nk from TMCOND'; STOP
         END IF
         DO jj = 1,iCOMP
            IF (Mko(k,jj).LT.-1.0d-10) THEN
               PRINT*, 'Negative Mk from TMCOND'; STOP
            END IF
         END DO
      END DO

!     UPDATE GAS-PHASE CONCENTRATION:
      Gcf(j) = &
      Gcf(j) + (SUM(Mkf(:,j)) - SUM(Mko(:,j)))
      
      IF (Gcf(j).LT.0.0) THEN
         IF ((ABS(Gcf(j)/SUM(Mkf(:,j))).LT.1e-8).OR.(ABS(Gcf(j)).LT.1e-8)) THEN
            Gcf(j) = 0.0
         ELSE
            PRINT*,'Negative Gc after TMCOND'; STOP
         END IF
      END IF
      
!     SWAP ARRAYS:
      DO k = 1,iBINS
         Nkf(k) = Nko(k)
         DO jj = 1,iCOMP
            Mkf(k,jj) = Mko(k,jj)
         END DO
      END DO
      
!     TIME LOOP CONTROL:
!     ====================================================================================
      CTR = CTR + 1; IF (CTR.LT.nSTEP) GOTO 100
      
      END DO LOOP1
      
      DO k = 1,iBINS
         Nk(k) = Nkf(k)
         DO jj = 1,iCOMP
            Mk(k,jj) = Mkf(k,jj)
         END DO
      END DO

      DO jj = 1,iCOMP-1
         Gc(jj) = Gcf(jj)
      END DO
      
      RETURN
      END SUBROUTINE






