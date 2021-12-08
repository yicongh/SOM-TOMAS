!     ===================================================================================
!             THIS SUBROUTINE INITIALIZES THE SOM GAS-PHASE ARRAYS FOR THE MODEL
!     ===================================================================================
      
      SUBROUTINE INIT_GPHASE(CLASS_NAMES, CLASS_DLVPS, &
                             PRECS_NAMES, PRECS_IPPMS, &
                             Temp, Pres, BoxVol)
      
      USE hSOM; USE hTOMAS; IMPLICIT NONE
      
!     DECLARATIONS:
!     ===================================================================================
!     SOM GRID NAMES AND
!     DLVP VALUES:
      CHARACTER*16     CLASS_NAMES(nSOMCLASS)
      DOUBLE PRECISION CLASS_DLVPS(nSOMCLASS)

!     PRECURSOR NAMES AND
!     EMISSION AMOUNTS [ppm]:
      CHARACTER*16     PRECS_NAMES(nSOMPRECS)
      DOUBLE PRECISION PRECS_IPPMS(nSOMPRECS)

!     CHAMBER TEMPERATURE [K],
!     PRESSURE [Pa] AND
!     CHAMBER VOLUME [cm3]:
      DOUBLE PRECISION Temp
      DOUBLE PRECISION Pres
      DOUBLE PRECISION BoxVol
      
!     SOM SPECIES ID,
!     DLVP VALUES AND
!     MOLECULAR WEIGHT OF CARBON BACKBONE [g mol-1]:
      CHARACTER*6      SOM_IDS(iORG)
      DOUBLE PRECISION SOM_DLVPS(iORG)
      DOUBLE PRECISION SOM_MWBONE(iORG)
      
!     FILE NAMES:
      CHARACTER*120 NAME1,NAME2,NAME3,NAME4
    
!     STRING HOLDER:
      CHARACTER*16 KEY
      
!     INDEX AND
!     COUNTERS:
      INTEGER INDEX,INDEX1,INDEX2
      INTEGER i,j,k,p,q

!     OFFSET OF ACTIVE SPECIES NAMES:
      INTEGER OFFSET
      PARAMETER(OFFSET = 13)
      
!     GAS CONSTANT [J mol-1 K-1]:
      DOUBLE PRECISION R
      PARAMETER(R = 8.314)
      
!     REFERENCE TEMPERATURE [K]:
      DOUBLE PRECISION T_ref
      PARAMETER(T_ref = 298.0)
      
!     ASSIGN SOM SPECIES NAMES:
!     ===================================================================================
      DO i = 1,nSOMPRECS
         SOMNAMES(i) = PRECS_NAMES(i)
      END DO
      
      p = nSOMPRECS + 1
      
      DO i = 1,nSOMCLASS
      DO j = 1,GRID_CMAX(i)
      DO k = 1,MIN(2*j,GRID_OMAX(i))
         WRITE(SOMNAMES(p),"((A8),'_',(I0.2),'_',(I0.2))") CLASS_NAMES(i),j,k; p = p + 1
      END DO
      END DO
      END DO
      
!     ASSIGN TOMAS SPECIES NAMES:
!     ===================================================================================
!     INORGANIC SPECIES:
      TOMASNAMES(1)       = 'SO4'
      TOMASNAMES(iCOMP-2) = 'POA'
      TOMASNAMES(iCOMP-1) = 'NH4'
      TOMASNAMES(iCOMP)   = 'H2O'
      
!     MONOMER SPECIES NAMES:
      TOMASNAMES(xORG_1ST:xMONO_LST) = SOMNAMES(nSOMPRECS+1:)

!     DIMER SPECIES:
      TOMASNAMES(xDIME_1ST:xORG_LST-1) = TOMASNAMES(xORG_1ST:xMONO_LST)
      
      DO i = 1,iMONO
         WRITE(TOMASNAMES(xDIME_1ST+i-1)(5:8),'(A4)') 'DIME'
      END DO
      
!     CARBON AND OXYGEN NUM. OF THE MONOMERS:
!     ===================================================================================
      CNUMBER(:) = 1
      ONUMBER(:) = 0
      
      DO i = 1,iORG
         READ(TOMASNAMES(xORG_1ST+i-1)(10:11),'(I2)') CNUMBER(i)
         READ(TOMASNAMES(xORG_1ST+i-1)(13:14),'(I2)') ONUMBER(i)
      END DO
      
!     ADAPTATION FOR HOMS:
      DO i = 1,nSOMCLASS
         INDEX1 = xGRID_LST(i) - xGRID_1ST(1) + 1
         INDEX2 = xGRID_LST(i) - xGRID_1ST(1) + xORG_1ST
         
         CNUMBER(INDEX1) = CNUMBER(INDEX1)-1
         ONUMBER(INDEX1) = MAX(CNUMBER(INDEX1),GRID_OMAX(i)+1)
      END DO
      
      CNUMBER(iORG) = 30
      ONUMBER(iORG) = 0
      
!     CALCULATE O:C RATIO FOR EACH SPECIES:
!     ===================================================================================
      DO i = 1,iORG
         ORG_O2C(i) = DBLE(ONUMBER(i))/DBLE(CNUMBER(i))
      END DO
      
!     DLVP VALUE FOR EACH MONOMER:
!     ===================================================================================      
      SOM_IDS(:) = TOMASNAMES(xORG_1ST:xORG_LST)(1:4); SOM_DLVPS(:) = 0.0
      
      DO i = 1,nSOMCLASS
         WHERE(SOM_IDS==CLASS_NAMES(i)(1:4))
            SOM_DLVPS = CLASS_DLVPS(i)
         END WHERE
      END DO
      
!     MOLECULAR WEIGHT OF THE PECIES:
!     ===================================================================================      
      ORG_MW = CNUMBER*12.0107 + ONUMBER*15.999 + (CNUMBER*2 + 2 - ONUMBER)*1.00794
      
!     GAS-PHASE CONC. FOR TOMAS AND SAPRC ARRAYS:
!     ===================================================================================
      SOMGC(:) = 0.0; GC(:) = 0.0
      
!     INITIAL CONC. IN SAPRC:      
      DO i = 1,nSOMPRECS
         WHERE(SOMNAMES==PRECS_NAMES(i))
            SOMGC = PRECS_IPPMS(i)
         END WHERE
      END DO
      
!     INITAL CONC. IN TOMAS:
      DO i = 1,iMONO
         GC(xORG_1ST+i-1) = &
         DBLE(SOMGC(nSOMPRECS+i))*1.0E-6*PRES/R/TEMP*ORG_MW(i)*1.0E-3*BOXVOL*1.0e-6
      END DO
      
!     THERMODYNAMIC PROPERTIES:
!     ===================================================================================
!     MOLECULAR WEIGHT OF CARBON BACKBONE:
      SOM_MWBONE = CNUMBER*12.0107 + (CNUMBER*2 + 2)*1.00794

!     SATURATION CONC.:
      Csat_Ref = 10.0**(-0.0337*SOM_MWBONE + 11.56 - ONUMBER*SOM_DLVPS)

!     FOR THE POA SPECIES:
      Csat_Ref(iORG) = 1e-6      
      Csat_Ref(iMONO-3) = 1e-6 !!!!!! ADDED HERE
      
!     SATURATION PRESSURE:
      Psat_Ref = Csat_Ref*R*T_ref/ORG_MW/1.0e6

!     ENTHALPY OF EVAPORATION:
      Hvap = -11.0*LOG10(Csat_Ref) + 131.0
      
!     SAVE DATA:
!     ===================================================================================            
      NAME1 = 'outputs/out.tomasnames'
      NAME2 = 'outputs/out.somnames'
      NAME3 = 'outputs/out.csat'
      NAME4 = 'outputs/out.o2c'

      OPEN(UNIT=345,FILE=NAME1,STATUS='unknown')
      WRITE(345,*) TOMASNAMES
      CLOSE(345)
      
      OPEN(UNIT=345,FILE=NAME2,STATUS='unknown')
      WRITE(345,*) SOMNAMES
      CLOSE(345)

      OPEN(UNIT=345,FILE=NAME3,STATUS='unknown')
      WRITE(345,*) CSAT_REF
      CLOSE(345)

      OPEN(UNIT=345,FILE=NAME4,STATUS='unknown')
      WRITE(345,*) ORG_O2C
      CLOSE(345)
      
!     SET TRIVIAL VALUES:
!     ===================================================================================
!     SO2 INITIAL CONC.:
      SO2 = 0.001*1e-3
      
      DO i = 1,iCOMP-1
         IF (GC(i).LT.1d-60) THEN
            GC(i) = 1d-60
         END IF
      END DO
      
      DO i = 1,nSOM
         IF (SOMGC(i).LT.1d-60) THEN
            SOMGC(i) = 1d-60
         END IF
      END DO
      
 300  RETURN
      END SUBROUTINE
