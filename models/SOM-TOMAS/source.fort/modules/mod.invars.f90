!     ==========================================================================
!                     THIS IS THE MODULE FOR THE INPUT VARIABLES
!     ==========================================================================

      MODULE hINVARS
            
      IMPLICIT NONE

!     SIMULATION NAME:
      CHARACTER(LEN=120),SAVE :: SIMUNAME
      
!     MODEL TIMESTEP [s] AND
!     OUTPUT FREQUENCY:
      REAL(8),SAVE :: mdt
      INTEGER,SAVE :: freq
      
!     TOTAL SIMULATION TIME [s]:
      REAL(8),SAVE :: t_end
      
!     SWITCHES [0 or 1]:
      INTEGER,SAVE :: mCOAG
      INTEGER,SAVE :: mVWL
      INTEGER,SAVE :: mPWL
      INTEGER,SAVE :: mHET
      INTEGER,SAVE :: mOLIG
      INTEGER,SAVE :: mENDB
      INTEGER,SAVE :: mGEN1
      INTEGER,SAVE :: mV2PWL
      INTEGER,SAVE :: mSORB
      INTEGER,SAVE :: mNUCL
      INTEGER,SAVE :: mWSORB
      
!     REACTOR CONDITIONS - PRESSURE [Pa],
!     TEMPERATURE [K],
!     RELATIVE HUMIDITY AND
!     REACTOR VOLUME [cm3]:
      REAL(8),SAVE :: Pres
      REAL(8),SAVE :: Temp
      REAL(8),SAVE :: RH
      REAL(8),SAVE :: BoxVol
      
!     VAPOR WALL LOSS RATE [s-1]:
      REAL(8),SAVE :: kvap_on
      
!     AEROSOL PROPERTIES - ACCOMMODATION COEFFICIENT,
!     SURFACE TENSION [N m-1],
!     BULK DIFFUSIVITY [m2 s-1],
!     CHEMICAL DECAY RATE [s-1] AND
!     NUCLEATION RATE [s-1]:
      REAL(8),SAVE :: Alpha
      REAL(8),SAVE :: ST
      REAL(8),SAVE :: Db
      REAL(8),SAVE :: kc
      REAL(8),SAVE :: NucRate
            
      REAL(8),SAVE :: OrgFrac
      REAL(8),SAVE :: RhoSeed

!     SOM GRID AND
!     PRECURSOR NAMES:
      CHARACTER(LEN=16),ALLOCATABLE,SAVE :: CLASS_NAMES(:)
      CHARACTER(LEN=16),ALLOCATABLE,SAVE :: PRECS_NAMES(:)
      
!     GRID DLVPS:
      REAL(8),ALLOCATABLE,SAVE :: CLASS_DLVPS(:)
            
!     PRECURSOR EMISSIONS:
      REAL(8),ALLOCATABLE,SAVE :: PRECS_IPPMS(:)

!     SOM PARAMETERS:
      REAL(8),ALLOCATABLE,SAVE :: SOM_PARAMS(:)

!     PRECURSOR kOH AND      
!     HOMS YIELD:
      REAL(8),ALLOCATABLE,SAVE :: PRECS_kOH(:)
      REAL(8),ALLOCATABLE,SAVE :: PRECS_fHOM(:)
      
!     PRECURSOR CARBON AND 
!     OXYGEN NUMBERS:
      INTEGER,ALLOCATABLE,SAVE :: PRECS_CNO(:)
      INTEGER,ALLOCATABLE,SAVE :: PRECS_ONO(:)
      
!     OLIGOMER FORMATION AND
!     DISSOCIATION RATES:
      REAL(8),SAVE ::  k_f
      REAL(8),SAVE ::  k_r

!     NUCLEATION PARAMETERIZATIONS:
      REAL(8),SAVE :: NUCL_A
      REAL(8),SAVE :: NUCL_B
      REAL(8),SAVE :: NUCL_C

!     BIN MASS BOUNDS:
      REAL(8),SAVE :: bMASS_MIN
      REAL(8),SAVE :: bMASS_MAX

!     PARTICLE WALL lOSS RATES:
      REAL(8),ALLOCATABLE,SAVE :: kpar(:)

!     OH TIMESERIES:
      REAL(8),ALLOCATABLE,SAVE :: i_OH_TS(:)

!     NUCL. RATE TIMESERIES:
      REAL(8),ALLOCATABLE,SAVE :: i_JNUC_TS(:)

!     OH UPTAKE COEFFICIENT:
      REAL(8),SAVE :: Gamma_OH
      
      END MODULE hINVARS
