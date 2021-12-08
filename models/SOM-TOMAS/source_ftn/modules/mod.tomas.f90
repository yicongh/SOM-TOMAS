!     ==========================================================================
!               THIS IS THE MODULE FOR THE TOMAS VARIABLES AND ARRAYS
!     ==========================================================================

      MODULE hTOMAS
      
      IMPLICIT NONE
      
!     NUMBER OF BINS,
!     COMPONENTS,
!     DIAGNOSTIC SPECIES
!     ORGANIC SPECIES AND
!     MONOMER SPECIES:
      INTEGER iBINS
      INTEGER iCOMP
      INTEGER iDIAG
      INTEGER iORG
      INTEGER iMONO
!FLAG1
      PARAMETER(iBINS=36, iCOMP=140, iDIAG=2, iORG=137)
      PARAMETER(iMONO=(iORG-1)/2)
      
!     INDICES OF SO4, NH4 AND H2O:
      INTEGER,PARAMETER :: xSO4 = 1
      INTEGER,PARAMETER :: xNH4 = iCOMP - 1
      INTEGER,PARAMETER :: xH2O = iCOMP
      
!     FIRST AND LAST INDICES OF ORGANICS:
      INTEGER,PARAMETER :: xORG_1ST = 2
      INTEGER,PARAMETER :: xORG_LST = xORG_1ST + iORG - 1
      
!     INDICES OF LAST MONOMER SPEICES AND
!     FIRST DIMER SPECIES:
      INTEGER,PARAMETER :: xMONO_LST = xORG_1ST + iMONO - 1
      INTEGER,PARAMETER :: xDIME_1ST = xMONO_LST + 1 
      
!     AEROSOL NUMBER [box-1],
!     MASS [kg box-1] AND
!     GAS-PHASE MASS [kg box-1]:
      REAL(8),SAVE :: Nk(iBINS)
      REAL(8),SAVE :: Mk(iBINS,iCOMP)
      REAL(8),SAVE :: Gc(iCOMP-1)

!     WALL-DEPOSITED AEROSOl NUMBER [box-1],
!     MASS [kg box-1] AND
!     GAS-PHASE MASS [kg box-1]:
      REAL(8),SAVE :: Nk_WALL(iBINS)
      REAL(8),SAVE :: Mk_WALL(iBINS,iCOMP)
      REAL(8),SAVE :: Gc_WALL(iCOMP-1)
      
!     MASS BOUNDS OF BINS [kg]:
      REAL(8),SAVE :: Xk(iBINS+1)
      
!     ENTHALPY OF EVAPORATION [J mol-1],
!     SATURATION CONCENTRATION [ug m-3] AND
!     SATURATION VAPOR PRESSURE [Pa]:
      REAL(8),SAVE :: Hvap(iORG)
      REAL(8),SAVE :: Csat_Ref(iORG)
      REAL(8),SAVE :: Psat_Ref(iORG)
      
!     CARBON AND
!     OXYGEN NUMBERS OF ORGANIC SPECIES:
      INTEGER,SAVE :: CNUMBER(iORG)
      INTEGER,SAVE :: ONUMBER(iORG)

!     MOLECULAR WEIGHT [g mol-1] AND 
!     O:C RATIO OF ORGANIC SPECIES:
      REAL(8),SAVE :: ORG_MW(iORG)
      REAL(8),SAVE :: ORG_O2C(iORG)
      
!     TOMAS SPECIES NAMES:
      CHARACTER(LEN=15),SAVE :: TOMASNAMES(iCOMP)
      
!     PARTICLE WALL LOSS RATES [s-1]:
      REAL(8),SAVE :: BETAS_PWL(iBINS)

!     HET. CHEM. AND OLIGOMERIZATIOIN
!     CONTRIBUTION TO PARTICLE-PHASE DECAY RATE [s-1]:
      REAL(8),SAVE :: kc_HET(iBINS,iCOMP) = 0.0
      REAL(8),SAVE :: kc_OLIG(iBINS,iCOMP) = 0.0
      
      END MODULE hTOMAS
