!     ==============================================================================
!               THIS SUBROUTINE STEPS FORWARD FOR HETEROGENEOUS CHEMISTRY
!     ==============================================================================
      
      SUBROUTINE STEP_HETCHEM(SOM_PARAMS,Gamma_OH,BoxVol,OH_in,Temp,Pres,mHET,dt)

      USE hTOMAS; USE hSOM; IMPLICIT NONE

!     LENGTH OF ARRAY:
      INTEGER LENGTH
      PARAMETER(LENGTH=nSOMCLASS*5)
      
!     SOM PARAMETERS:
      DOUBLE PRECISION SOM_PARAMS(LENGTH)
      
      DOUBLE PRECISION Temp
      DOUBLE PRECISION Pres
      INTEGER mHET

!     SOM PARAMETERS,
!     CHAMBER VOLUME [cm3],
!     OH CONC. [ppm] AND
!     TIMESTEP [s]:
      DOUBLE PRECISION PARAMS_IN(5)
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION OH_in
      DOUBLE PRECISION dt
      
      INTEGER nCELLS
      INTEGER INDEX_1ST

      INTEGER CMAX_IN
      INTEGER OMAX_IN

!     ARRAY FOR AEROSOL-PHASE SOM SPECIES [mole]:
      DOUBLE PRECISION,ALLOCATABLE :: AECONC(:)

!     MASS OF ORGANIC AND
!     SEED SPECIES IN EACH BIN [kg]:
      DOUBLE PRECISION mass_org
      DOUBLE PRECISION mass_ams

!     DENSITY OF ORGANIC AND
!     SEED MATERIAL [kg m-3]:
      DOUBLE PRECISION rho_org
      DOUBLE PRECISION rho_ams
      PARAMETER(rho_org = 1400.)
      PARAMETER(rho_ams = 1770.)

!     VOLUME [m3] AND
!     DIAMETER [m] OF SINGLE PARTICLES:
      DOUBLE PRECISION pp_volm
      DOUBLE PRECISION pp_diam
      
!     AVE. MW IN BIN [g mol-1]:
      DOUBLE PRECISION MW_AVE

!     OH IN [mole]:
      DOUBLE PRECISION OH

!     EFFECTIVE KOH FOR HET. CHEM. [mole-1 s-1]:
      DOUBLE PRECISION kOH_EFF

!     KNUDSEN NUMBER AND
!     FUCHS-SUTUGIN FACTOR:
      DOUBLE PRECISION Kn
      DOUBLE PRECISION FS

!     OH DIFFUSIVITY [m2 s-1],
!     OH UPTAKE COEFFICIENT,
!     OH ROOT MEAN VELOCITY [m s-1] AND
!     OH MEAN FREE PATH [m]:
      DOUBLE PRECISION D_OH
      PARAMETER(D_OH = 1.38e-5*44./17.)
      
      DOUBLE PRECISION Gamma_OH
!      PARAMETER(Gamma_OH = 1.0)
      
      DOUBLE PRECISION RMS_OH
      PARAMETER(RMS_OH = SQRT(3.*8.314*298./(17./1000.)))
      
      DOUBLE PRECISION MFP_OH
      PARAMETER(MFP_OH = 3.0*D_OH/RMS_OH)

!     PARTICLE SURFACE ACCOMMODATION COEFF.:
      DOUBLE PRECISION Alpha
      PARAMETER(Alpha = 1.0)

!     PI AND
!     AVOGADRO NUMBER [# mole-1]:
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.1415926)
      
      DOUBLE PRECISION NA
      PARAMETER(NA = 6.022e23)

!     COUNTERS:
      INTEGER i,j,k,jj,kk,n

!     DUMMY VARIABLE:
      DOUBLE PRECISION void

!     HET. CHEM. MODULE:
      DOUBLE PRECISION,EXTERNAL :: het_chem

      DOUBLE PRECISION Mk_NEXT(iBINS,iCOMP)

      DOUBLE PRECISION Mk_f(iBINS,iCOMP)      
      DOUBLE PRECISION Mk_o(iBINS,iCOMP)

      DOUBLE PRECISION Nk_f(iBINS)
      DOUBLE PRECISION Nk_o(iBINS)      
      
      DOUBLE PRECISION dMASS(iBINS)

      DOUBLE PRECISION X0(iBINS)
      DOUBLE PRECISION X1(iBINS)

      DOUBLE PRECISION TAU(iBINS)

      INTEGER MASK(iBINS)

      INTEGER nSTEPS
      PARAMETER(nSTEPS = 1)

      DOUBLE PRECISION mmdt
      DOUBLE PRECISION bMASS
      DOUBLE PRECISION FRAC
      
!     SKIPPING CRITERIA:
!     ==============================================================================
!     IF HET. CHEM. IS OFF:
      IF (mHET.EQ.0) GOTO 900
      
!     IF ORGANIC MASS IS SMALL:
      IF (SUM(Mk(:,xORG_1ST:xORG_LST))/BoxVol*1e6*1e9.LT.0.1) THEN ! 0.1
         GOTO 900
      END IF
      
!     STEP FOR HETEROGENEOUS CHEMISTRY:
!     ==============================================================================      
!     INTERNAL TIMESTEP:
      mmdt = dt/FLOAT(nSTEPS)
      
!     LOOP OVER INTERNAL STEPS:
      LOOP1: DO n = 1,nSTEPS

      DO k = 1,iBINS
         DO j = 1,iCOMP
            Mk_NEXT(k,j) = Mk(k,j)
         END DO
      END DO
      
!     LOOP OVER BINS:
      DO k = 1,iBINS
         IF (Nk(k).LT.100.0) GOTO 101

!        CALCULATE THE AVE. SIZE OF THE BIN:
         mass_org = SUM(Mk(k,xORG_1ST:xORG_LST))
         mass_ams = Mk(k,xSO4)

         pp_volm = ((mass_org/rho_org) + (mass_ams/rho_ams))/Nk(k)
         pp_diam = (6.*pp_volm/PI)**(1./3.)

!         IF (pp_diam*1e9.LT.10.0) GOTO 101
         IF ((mass_org/(mass_org + mass_ams)).LT.1e-2) GOTO 101
         
!        CALCULATE THE KNUDSEN NUMBER:
         Kn = 2.0*MFP_OH/pp_diam

!        CALCULATE THE FUCHS-SUTUGIN FACTOR:
         FS = (1 + Kn)/(1 + Kn + 0.75/Kn)

!        CALCULATE AVE. MW:
         MW_AVE = 0.0
         
         DO j = 1,iORG
            MW_AVE = MW_AVE + ORG_MW(j)*Mk(k,xORG_1ST+j-1)/mass_org
         END DO
         
!        CALCULATE THE EFFECTIVE KOH [mole-1 s-1]:
         kOH_EFF = Gamma_OH*PI*(pp_diam**2.)*(Nk(k)/(BoxVol*1e-6))* &
                   FS*RMS_OH/4./(mass_org/(MW_AVE/1000.))
         
!        CONVERT OH TO [mole]:
         OH = OH_in*(Pres/8.314/Temp)*1e-6*(BoxVol*1e-6)
         
!        REGULATE REACTION RATE:
         IF (kOH_EFF*OH*mmdt.GT.3.0) THEN
            OH = 3.0/kOH_EFF/mmdt
         END IF
         
!        LOOP OVER SOM CLASSES:
         jj = 1
         kk = 1

         DO i = 1,nSOMCLASS
!           SELECT SOM PARAMETERS:
            PARAMS_IN = SOM_PARAMS(5*i-4:5*i)
         
!           NUMBER OF CELLS:
            nCELLS = uSOMCELLS(i)         
            
!           READ MASS CONCENTRATION:
            ALLOCATE(AECONC(nCELLS))
            
            INDEX_1ST = xGRID_1ST(i)-xGRID_1ST(1)+xORG_1ST

            DO j = 1,nCELLS
               AECONC(j) = Mk(k,INDEX_1ST+j-1)/(ORG_MW(jj)/1000.)*1e20
               jj = jj + 1
            END DO
            
!           MAX CARBON AND
!           OXYGEN NUMBERS:
            CMAX_IN = GRID_CMAX(i)
            OMAX_IN = GRID_OMAX(i)
            
!           STEP FOR HETEROGENEOUS CHEM.:
            void = het_chem(AECONC,PARAMS_IN,OH,kOH_EFF, &
                            Temp,Pres,mmdt, &
                            nCELLS,CMAX_IN,OMAX_IN)
            
            DO j = 1,nCELLS
               Mk_NEXT(k,INDEX_1ST+j-1) = &
               AECONC(j)*(ORG_MW(kk)/1000.)*1e-20
               kk = kk + 1
            END DO
            
            DEALLOCATE(AECONC)

         END DO
 101     CONTINUE
      END DO

!     CALCULATE HET. CHEM. CONTRIBUTION TO PARTICLE-PHASE DECAY:
!     ==============================================================================
      IF (n.EQ.nSTEPS) THEN
!        REFRESH THE ARRAY:
         kc_HET(:,:) = 0.0

!        PSEUDO-FIRST-ORDER DECAY:         
         DO kk = 1,iBINS
            DO jj = 1,iCOMP
               kc_HET(kk,jj) = LOG(MIN(Mk_NEXT(kk,jj)/Mk(kk,jj),1.0))/(-mmdt)
            END DO
         END DO
         
      END IF
      
!     REMAP NUMBER AND MASS:
!     ==============================================================================
      MASK(:) = 1

      DO k = 1,iBINS
         IF (Nk(k).LT.5.) MASK(k) = 0
      END DO
      
      X0 = SUM(Mk,2)/Nk
      X1 = SUM(Mk_NEXT,2)/Nk

      TAU = (3./2.)*(X1**(2./3.) - X0**(2./3.))

      DO i = 1,iBINS
         Nk_f(i) = Nk(i)
         DO k = 1,iCOMP
            Mk_f(i,k) = Mk(i,k)
         END DO
      END DO
      
      CALL TMCOND(TAU,Xk,Mk_f,Nk_f,Mk_o,Nk_o,xORG_LST+1,MASK)
      
!      DO k = 1,iBINS
!         bMASS = SUM(Mk_o(k,:))
!         DO j = 1,iCOMP-2
!            FRAC = Mk_NEXT(k,j)/SUM(Mk_NEXT(k,1:iCOMP-2))
!            Mk_o(k,j) = bMASS*FRAC
!         END DO
!         Mk_o(k,iCOMP-1) = 0.d0
!      END DO
      
      DO i = 1,iBINS
         Nk(i) = Nk_o(i)
         DO k = 1,iCOMP-2
            Mk(i,k) = Mk_NEXT(i,k)*Mk_o(i,k)/Mk_f(i,k)  !Mk_o(i,k)  !!! IMPORTANT HERE !!!!!!!!!!!!!!!!!!!!
         END DO
      END DO

      END DO LOOP1
      
 900  RETURN
      END SUBROUTINE
