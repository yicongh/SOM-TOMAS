C ==============================================================================
C This function calculates the aerosol bulk diffusivity based on its composition
C (monomers, dimers and water).
C ==============================================================================

      DOUBLE PRECISION FUNCTION CALC_DB(Temp,Pres,RH,BoxVol,timer)

      USE hTOMAS; USE hSOM; IMPLICIT NONE
      
C      INCLUDE "../headers/h_tomas_arrays.f"
C      INCLUDE "../headers/h_som_arrays.f"

      DOUBLE PRECISION ofrac
      DOUBLE PRECISION Temp            ! [K]  Temperature
      DOUBLE PRECISION Pres            ! [Pa] Pressure
      DOUBLE PRECISION RH              !      Relative humidity
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION timer
      
      DOUBLE PRECISION Mass_Org          ! [kg] Total mass of organics (monomers + dimers)
      DOUBLE PRECISION Frac_Orgs(iMONO+1)!(iORG) !      Fraction of each organic species
      DOUBLE PRECISION Frac_Orgs2(iORG)
      
      DOUBLE PRECISION Tg_Orgs(iMONO+1) !(iORG)   ! [K]  Glass-transition temperature of each organic species
      DOUBLE PRECISION Tg_Org           ! [K]  Ave. glass-transition temperature of total organics
      DOUBLE PRECISION Tg               ! [K]  Ave. glass-transition temperature of organics and water

      DOUBLE PRECISION Mass_W          ! [kg] Mass of water
      DOUBLE PRECISION Frac_W          !      Fraction of water

      DOUBLE PRECISION eta             ! [Pa s] Aerosol viscosity at current temperature
      DOUBLE PRECISION eta_0           ! [Pa s] Aerosol viscosity at infinite temperature
      DOUBLE PRECISION eta_c           ! [Pa s] Crossover viscosity (Evoy et al., 2019)
      PARAMETER(eta_0 = 1e-5)
      PARAMETER(eta_c = 1e-3)         

      DOUBLE PRECISION zeta            !        Correction factor (Evoy et al., 2019)
      PARAMETER(zeta = 0.93)

      DOUBLE PRECISION MW_ave          ! [g mol-1] Ave. MW of the condensing species weighted by mass
      DOUBLE PRECISION d_mo            ! [m]       Ave. diameter of condensing molecules
      DOUBLE PRECISION Cc              !           Slip-correction factor

      DOUBLE PRECISION MW_ave2
      
      DOUBLE PRECISION PI
      PARAMETER(PI = 3.1415926)

      DOUBLE PRECISION Na              ! [mol-1] Avogadro constant
      PARAMETER(Na = 6.022e23)
      
      DOUBLE PRECISION kB              ! [m2 kg s-2 K-1] Boltzmann constant
      PARAMETER(kB = 1.38e-23)    

      DOUBLE PRECISION k_GT            ! Gordon-Taylor constant (Shiraiwa, 2017)
      PARAMETER(k_GT = 2.5)

      DOUBLE PRECISION Tg_w            ! [K] Glass-transition temperature of water
      PARAMETER(Tg_w = 136.0)
      
      DOUBLE PRECISION D_frag          ! Fragility of aerosol liquid
      PARAMETER(D_frag = 10.0)

      DOUBLE PRECISION rho_org
      PARAMETER(rho_org = 1400.0)      ! [kg m-3] Density of organic aerosol

      DOUBLE PRECISION rho_w           ! [kg m-3] Density of water
      PARAMETER(rho_w = 1000.0)
      
      DOUBLE PRECISION beta            ! Ratio of Tg to To (Angell, 1995)
      PARAMETER(beta = 0.7966)

      DOUBLE PRECISION mfp             ! [m] Mean free path of air molecules
      PARAMETER(mfp = 66e-9)           

      DOUBLE PRECISION kappa           ! Kappa value of the organics
      PARAMETER(kappa = 0.1)
      
      DOUBLE PRECISION A,B,C,D,E       ! Parameterizations of glass-transition temperature model (Shiraiwa et al., 2017)
      PARAMETER(A = -21.57)
      PARAMETER(B = 1.51)
      PARAMETER(C = -1.7E-3)
      PARAMETER(D = 131.4)
      PARAMETER(E = -0.25)
      
      INTEGER i,j                      ! Loop counters

      DOUBLE PRECISION Mass_HOMS
      DOUBLE PRECISION Mass_LESS

      COMMON /TG_OUT/ Tg

      DOUBLE PRECISION MW_DIME
      DOUBLE PRECISION O2C_DIME

      INTEGER INDEX
      
C     CALCULATE AVE. GLASS-TRANSITION TEMPERATURE OF THE ORGANICS:
C     ==========================================================================
C     GLASS-TRANSITION TEMPERATURE OF EACH SPECIES:
      DO i = 1,iMONO
         Tg_Orgs(i) = A + B*ORG_MW(i) + C*(ORG_MW(i)**2.) +D*ORG_O2C(i)+ 
     &                E*ORG_MW(i)*ORG_O2C(i)
      END DO
      
C     AVE. MOLECULAR WEIGHT OF DIMERS:
      MW_DIME = 0.0
      
      DO i = 1,iMONO
         MW_DIME = 
     &   MW_DIME + ORG_MW(i)*SUM(Mk(:,xDIME_1ST+i-1))/
     &             SUM(Mk(:,xDIME_1ST:xORG_LST))*2.
      END DO
      
C     AVE. O:C RATIO OF DIMERS:
      O2C_DIME = 0.0

      DO i = 1,iMONO
         O2C_DIME = 
     &   O2C_DIME + ORG_O2C(i)*SUM(Mk(:,xDIME_1ST+i-1))/
     &              SUM(Mk(:,xDIME_1ST:xORG_LST))
      END DO
      
      Tg_Orgs(iMONO+1) = A + B*MW_DIME + C*(MW_DIME**2.) + D*O2C_DIME + 
     &                   E*MW_DIME*O2C_DIME
      
C     MASS FRACTION OF EACH SPECIES:
      Mass_ORG = SUM(Mk(:,xORG_1ST:xORG_LST))
      
      DO i = 1,nSOMCLASS   
         INDEX = xGRID_LST(i) - xGRID_1ST(1) + 1
         
C         Mass_ORG = 
C     &   Mass_ORG - SUM(Mk(:,xORG_1ST+INDEX-1))
      END DO
      
      DO i = 1,iMONO
         Frac_Orgs(i) = SUM(Mk(:,xORG_1ST+i-1))/Mass_ORG
      END DO

      Frac_Orgs(iMONO+1) = SUM(Mk(:,xMONO_LST+1:xORG_LST))/Mass_ORG
      
C     SKIP MASS FRACTION FOR HOMS:
      DO i = 1,nSOMCLASS
         
         INDEX = xGRID_LST(i) - xGRID_1ST(1) + 1
         
C         Frac_Orgs(INDEX) = 0.0
         
      END DO
      
C     AVE. GLASS-TRANSITION TEMPERATURE:
      Tg_Org = 0.0
      
      DO i = 1,iMONO+1
         Tg_Org = 
     &   Tg_Org + Tg_Orgs(i)*Frac_Orgs(i)
      END DO
      
C     ACCOUNT FOR WATER CONTENT IN THE AEROSOL PHASE:
C     ==========================================================================      
C     MASS FRACTION OF WATER:      
      Mass_W = Mass_Org*(rho_w/rho_org)*(kappa/(1.0/RH - 1.0))
      Frac_W = Mass_W/(Mass_W + Mass_Org)
      
C     OVERALL EFFECTIVE GLASS-TRANSITION TEMPERATURE:
      Tg = (Frac_W*Tg_W + (1./K_GT)*(1. - Frac_W)*Tg_Org) /
     &     (Frac_W + (1./K_GT)*(1. - Frac_W))

C     CALCULATE BULK DIFFUSIVITY:
C     ==========================================================================
C     AVE. MOLECULAR WEIGHT OF THE CONDENSING SPECIES:
      MW_ave = 0.0
      
      DO i = 1,iMONO
         MW_ave = 
     &   MW_ave + ORG_MW(i)*Frac_Orgs(i)/(1.-Frac_Orgs(iMONO+1))
      END DO
       
C     CALCULATE AEROSOL VISCOSITY AT CURRENT TEMPERATURE:
      eta = eta_0*EXP(D_frag*beta/((Temp/Tg) - beta))

      IF (eta.GT.1e12) eta = 1e12
      
C     AVE. DIAMETER OF THE CONDENSING MOLECULES:
      d_mo = ((6./PI)*(MW_ave/1000.)/(rho_org*Na))**(1./3.)
      
C     BULK DIFFUSIVITY:
      CALC_DB = kB*Temp/(6.*PI*d_mo*eta_c)*(eta_c/eta)**(zeta)
      
 100  RETURN
      END FUNCTION
