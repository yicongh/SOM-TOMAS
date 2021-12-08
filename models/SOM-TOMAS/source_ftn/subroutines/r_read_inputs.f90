!     ===================================================================================
!                   THIS SUBROUTINE READS THE VALUES OF THE INPUT VARIABLES
!     ===================================================================================
      
      SUBROUTINE READ_INPUTS(nSOMCLASS, nSOMPRECS, iBINS)
      
      USE hINVARS; IMPLICIT NONE
      
!     DECLARATIONS:
!     ===================================================================================
!     NUMBER OF SOM CLASSES,
!     PRECURSORS, BINS:
      INTEGER nSOMCLASS
      INTEGER nSOMPRECS
      INTEGER iBINS

!     NUMBER OF TIMESTEPS:
      INTEGER nSTEPS
      
!     ALLOCATE ARRAYS:
!     ===================================================================================
      ALLOCATE(CLASS_NAMES(nSOMCLASS))
      ALLOCATE(CLASS_DLVPS(nSOMCLASS))

      ALLOCATE(PRECS_NAMES(nSOMPRECS))
      ALLOCATE(PRECS_IPPMS(nSOMPRECS))
      
      ALLOCATE(PRECS_kOH(nSOMPRECS))
      ALLOCATE(PRECS_CNO(nSOMPRECS))
      ALLOCATE(PRECS_ONO(nSOMPRECS))
      ALLOCATE(PRECS_fHOM(nSOMPRECS))
      
      ALLOCATE(SOM_PARAMS(nSOMCLASS*5))
      
      ALLOCATE(kpar(iBINS))

!     READ THE INPUT VARIABLES:
!     ===================================================================================
!     OPEN THE INPUT FILE: 
      OPEN(UNIT=30,FILE='inputs/in.vars',STATUS='OLD')
      
      READ(30,*) SIMUNAME

      READ(30,*) mdt
      READ(30,*) freq
      READ(30,*) t_end
      
      READ(30,*) mCOAG
      READ(30,*) mVWL
      READ(30,*) mPWL
      READ(30,*) mHET
      READ(30,*) mENDB
      READ(30,*) mGEN1
      READ(30,*) mV2PWL
      READ(30,*) mSORB
      READ(30,*) mNUCL
      READ(30,*) mWSORB
      READ(30,*) mOLIG
      
      READ(30,*) Pres
      READ(30,*) Temp
      READ(30,*) RH
      READ(30,*) BoxVol
      READ(30,*) kvap_on
      
      READ(30,*) Alpha
      READ(30,*) Db

      READ(30,*) k_f
      READ(30,*) k_r

      READ(30,*) NUCL_A
      READ(30,*) NUCL_B
      READ(30,*) NUCL_C

      READ(30,*) CLASS_NAMES
      READ(30,*) CLASS_DLVPS

      READ(30,*) PRECS_NAMES
      READ(30,*) PRECS_IPPMS
      
      READ(30,*) SOM_PARAMS
      READ(30,*) PRECS_fHOM
      READ(30,*) PRECS_kOH
      
      READ(30,*) PRECS_CNO
      READ(30,*) PRECS_ONO

      READ(30,*) bMASS_MIN
      READ(30,*) bMASS_MAX
      
      READ(30,*) Gamma_OH

      CLOSE(30)

!     READ THE PARTICLE WALL LOSS RATES:
!     ===================================================================================
      OPEN(UNIT=30,FILE='inputs/in.PWL',STATUS='OLD')
      READ(30,*) kpar
      CLOSE(30)
      
!     READ THE OH PROFILE:
!     ===================================================================================
      nSTEPS = INT(t_end/mdt); ALLOCATE(i_OH_TS(nSTEPS))

      OPEN(UNIT=30,FILE='inputs/in.OH_PROFILE',STATUS='OLD')
      READ(30,*) i_OH_TS
      CLOSE(30)

!     READ THE OH PROFILE:
!     ===================================================================================
      ALLOCATE(i_JNUC_TS(nSTEPS))

      OPEN(UNIT=30,FILE='inputs/in.JNUC',STATUS='OLD')
      READ(30,*) i_JNUC_TS
      CLOSE(30)

      RETURN
      END SUBROUTINE
