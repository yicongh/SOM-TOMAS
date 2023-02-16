!     ==============================================================================
!                THIS SUBROUTINE STEPS FORWARD FOR GAS-PHASE CHEMISTRY
!     ==============================================================================

      SUBROUTINE STEP_GASCHEM(SOM_PARAMS,PRECS_kOH,PRECS_CNO,PRECS_ONO, &
                              Pres,Temp,BoxVol,OH,timer,G1ONLY,dt)

      USE hTOMAS; USE hSOM; IMPLICIT NONE
      
!     DECLARATIONS:
!     ==============================================================================
!     LENGTH OF ARRAY:
      INTEGER LENGTH
      PARAMETER(LENGTH=nSOMCLASS*5)
      
!     SOM PARAMETERS:
      DOUBLE PRECISION SOM_PARAMS(LENGTH)

!     PRECURSOR kOH [cm3 s-1]:
      DOUBLE PRECISION PRECS_kOH(nSOMPRECS)

!     PRECURSOR CARBON AND
!     OXYGEN NUMBERS:
      INTEGER PRECS_CNO(nSOMPRECS)
      INTEGER PRECS_ONO(nSOMPRECS)

      DOUBLE PRECISION Pres
      DOUBLE PRECISION Temp
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION OH
      DOUBLE PRECISION timer
      DOUBLE PRECISION dt
      
!     SWITCH FOR 1ST-GEN. ONLY:
      INTEGER G1ONLY

      DOUBLE PRECISION R
      PARAMETER(R = 8.314)
      
      INTEGER i,j,k
      INTEGER INDEX
      DOUBLE PRECISION dummy
      
      DOUBLE PRECISION,EXTERNAL :: GAS_CHEM
      
      INTEGER nPRECS,nCELLS,i1_PRECS,i1_CELLS,CNO(2),ONO(2)
      DOUBLE PRECISION kOH(2)
      
      INTEGER CMAX_IN
      INTEGER OMAX_IN
      
      DOUBLE PRECISION PARAMS_IN(5)
      
      DOUBLE PRECISION,ALLOCATABLE :: kOH_IN(:)
      INTEGER,ALLOCATABLE :: CNO_IN(:)
      INTEGER,ALLOCATABLE :: ONO_IN(:)

      IF (timer.LT.0.0) GOTO 100
      
!     UPDATE GAS-PHASE ARRAY:
!     ==============================================================================
      INDEX = xGRID_1ST(1)

      DO i = 1,iMONO
         SOMGC(INDEX+i-1) = 1.0E6*GC(i+xORG_1ST-1)*1.0E3/ORG_MW(i) / &
              (BOXVOL/1.0E6) / &
              (PRES/R/TEMP)
      END DO

!     STEP FOR GAS-PHASE CHEMISTRY:
!     ==============================================================================
      k = 1

      DO i = 1,nSOMCLASS
!        SELECT SOM PARAMETERS:
         PARAMS_IN = SOM_PARAMS(5*i-4:5*i)

!        NUMBER OF PRECURSORS AND CELLS:
         nPRECS = uSOMPRECS(i)
         nCELLS = uSOMCELLS(i)

!        INDEX OF FIRST PRECURSOR AND CELL;
!        MINUS 1 TO GO FROM FORTRAN TO C++:
         i1_PRECS = xPREC_1ST(i) - 1
         i1_CELLS = xGRID_1ST(i) - 1

!        PRECURSOR kOH,
!        CARBON NUMBER AND
!        OXYGEN NUMBER:
         ALLOCATE(kOH_IN(nPRECS))
         ALLOCATE(CNO_IN(nPRECS))
         ALLOCATE(ONO_IN(nPRECS))

         DO j = 1,nPRECS
            kOH_IN(j) = PRECS_kOH(k)
            CNO_IN(j) = PRECS_CNO(k)
            ONO_IN(j) = PRECS_ONO(k)
            k = k + 1
         END DO

!        MAX CARBON AND
!        OXYGEN NUMBERS:
         CMAX_IN = GRID_CMAX(i)
         OMAX_IN = GRID_OMAX(i)

         dummy = GAS_CHEM(SOMGC,PARAMS_IN,OH,Temp,Pres,dt, &
                       nPRECS,nCELLS,i1_PRECS,i1_CELLS, &
                       CNO_IN,ONO_IN,kOH_IN, &
                       CMAX_IN,OMAX_IN,G1ONLY) 
         
         DEALLOCATE(kOH_IN)
         DEALLOCATE(CNO_IN)
         DEALLOCATE(ONO_IN)


      END DO

!     UPDATE GAS-PHASE ARRAY:
!     ==============================================================================
      INDEX = xGRID_1ST(1)
      
      DO i = 1,iMONO
         GC(i+xORG_1ST-1) = SOMGC(INDEX+i-1)* &
             (PRES/R/TEMP)*(BOXVOL/1.0E6) / &
             (1.0E6*1.0E3/ORG_MW(i))         
      END DO
      
 100  CONTINUE
      RETURN
      END SUBROUTINE
