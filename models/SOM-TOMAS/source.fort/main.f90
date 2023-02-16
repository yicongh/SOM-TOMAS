!     ===================================================================================
!                      THIS IS THE DRIVER PROGRAM FOR THE SOM-TOMAS MODEL
!     ===================================================================================
      
      PROGRAM BOX
        
      USE hMAIN; IMPLICIT NONE

!     INPUT AND OUTPUT DIRECTORIES:
!     ===================================================================================
      WRITE(*,*) "SPECIFY INPUT FOLDER: "
      READ(*,'(a)') PATHi

      WRITE(*,*) "SPECIFY OUTPUT FOLDER: "
      READ(*,'(a)') PATHo
      
!     READ INPUT VARIABLES:
!     ===================================================================================
      CALL READ_INPUTS(nSOMCLASS, nSOMPRECS, iBINS, PATHi)
                  
!     INITIALIZE GAS-PHASE ARRAYS:
!     ===================================================================================
      CALL INIT_GPHASE(CLASS_NAMES, CLASS_DLVPS, &
                       PRECS_NAMES, PRECS_IPPMS, &
                       Temp, Pres, BoxVol, PATHo)
      
!     INITIALIZE PARTICLE-PHASE ARRAYS:
!     ===================================================================================      
      CALL INIT_PPHASE(bMASS_MIN, bMASS_MAX, mSORB, PATHi)
            
!     INITIALIZE THE OUTPUT FILES:
!     ===================================================================================
      CALL INIT_OUTFILES(SIMUNAME, PATHo)
      
!     ITERATIONS START FROM HERE:
!     ===================================================================================      
      timer = 0.0; iCTR = 1
      
      CALL WRITE_OUTPUTS(timer,Db,OH,BoxVol,JNUC)
      
      LOOP0: DO WHILE (NINT(timer).LT.t_end)

      WRITE(*,"('time = ',1X,(F10.1))") timer
         
      LOOP1: DO i = 1,freq
      
!     STEP FOR NUCLEATION:
!     ===================================================================================
      IF (mNUCL.EQ.0) THEN
         JNUC = CALC_JNUC(NUCL_A,NUCL_B,NUCL_C,Pres,Temp,iCTR)
      ELSE
         JNUC = i_JNUC_TS(iCTR)
      END IF
      
      CALL STEP_NUCL(JNUC,BoxVol,timer,mdt)
      
!     STEP FOR PARTICLE AND VAPOR WALL LOSS:
!     ===================================================================================
      CALL STEP_WALLLOSS(kvap_on,kpar,mVWL,mPWL,mWSORB,Temp,BoxVol,timer,mdt)
      
!     STEP FOR COAGULATION:
!     ===================================================================================
      CALL STEP_COAG(Temp,Pres,BoxVol,mCOAG,mdt)
      
!     STEP FOR GAS-PHASE CHEMISTRY:
!     ===================================================================================
      OH = i_OH_TS(iCTR); OH = OH*1.0E6/Na/(Pres/R/Temp)*1.0E6

      CALL STEP_AUTOOX(PRECS_fHOM,PRECS_kOH,OH,Temp,Pres,BoxVol,timer,mdt)
      
      CALL STEP_GASCHEM(SOM_PARAMS, &
                        PRECS_kOH,PRECS_CNO,PRECS_ONO, &
                        Pres,Temp,BoxVol,OH,timer,mGEN1,mdt)
      
!     STEP FOR HETEROGENEOUS CHEMISTRY:
!     ===================================================================================
      CALL STEP_HETCHEM(SOM_PARAMS,Gamma_OH,BoxVol,OH,Temp,Pres,mHET,mdt)

!     STEP FOR AEROSOL-PHASE CHEMISTRY:
!     ===================================================================================
      CALL STEP_AEROCHEM(k_f,k_r,mOLIG,mdt,timer)

!     STEP FOR CONDENSATION:
!     ===================================================================================
!!      IF (ENDB.EQ.1) THEN
!!         Db = CALC_DB(Temp,Pres,RH,BoxVol,timer)
!!      END IF
      
      CALL STEP_SOACOND(Pres,Temp,RH,BoxVol,Alpha,Db,mdt)
      
      timer = &
      timer + mdt
      
      iCTR = &
      iCTR + 1
      
      END DO LOOP1
      
!     WRITE OUTPUTS:
!     ===================================================================================
      CALL WRITE_OUTPUTS(timer,Db,OH,BoxVol,JNUC)

      END DO LOOP0
      
      PRINT*, 'FINISHED'
      END PROGRAM

