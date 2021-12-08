!     ==========================================================================
!                THIS IS THE MODULE FOR THE SOM VARIABLES AND ARRAYS
!     ==========================================================================

      MODULE hSOM
      
      IMPLICIT NONE
      
!     NUMBER OF SOM CLASSES AND
!     PRECURSORS:
!FLAG1
      INTEGER,PARAMETER :: nSOMCLASS = 1
      INTEGER,PARAMETER :: nSOMPRECS = 1
      
!     NUMBER OF SOM PRECURSORS IN EACH CLASS:
!FLAG2
      INTEGER,PARAMETER :: uSOMPRECS(nSOMCLASS) = (/1/)
      INTEGER,PARAMETER :: uSOMCELLS(nSOMCLASS) = (/68/)
      
!     INDICES OF FIRST PREC. AND GRID SPECIES:
!FLAG3
      INTEGER,PARAMETER :: xPREC_1ST(nSOMCLASS) = (/1/)
      INTEGER,PARAMETER :: xGRID_1ST(nSOMCLASS) = (/2/)
      INTEGER,PARAMETER :: xGRID_LST(nSOMCLASS) = (/69/)

!     SIZE OF SOM ARRAYS:
!FLAG4
      INTEGER,PARAMETER :: nSOM = 69
      
!     GAS-PHASE CONC.:
      REAL(8),SAVE :: SOMGC(nSOM)
      
!     GAS-PHASE SPECIES NAMES:
      CHARACTER(LEN=20) :: SOMNAMES(nSOM)
      
!     NUMBER OF SOM PARAMETERS:
!FLAG5
      INTEGER,PARAMETER :: nSOMPARAM = 5
      
!     MAX CARBON AND OXYGEN NUMBERS:
!FLAG6
      INTEGER,PARAMETER :: GRID_CMAX(nSOMCLASS) = (/11/)
      INTEGER,PARAMETER :: GRID_OMAX(nSOMCLASS) = (/7/)

      REAL(8),SAVE :: SO2      
      
      END MODULE hSOM
