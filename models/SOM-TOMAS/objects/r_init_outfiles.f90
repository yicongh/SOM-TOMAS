!     ==================================================================
!              THIS SUBROUTINE INITIALIZES THE OUTPUT FILES
!     ==================================================================
      
      SUBROUTINE INIT_OUTFILES(SIMUNAME)
      
      USE hSOM; USE hTOMAS; USE hOFILES; IMPLICIT NONE

!     DECLARATIONS:
!     ==================================================================
!     SIMULATION NAME:      
      CHARACTER(120) SIMUNAME

!     NAME LENGTH:      
      INTEGER LENGTH
      
!     FUNCTION TO COUNT LENGTH:
      INTEGER,EXTERNAL :: GET_LEN

!     LOOP COUNTERS:
      INTEGER i,j,k

!     FORMAT STRING:
      CHARACTER(120) FMT

!     ASSIGN OUTFILE PATHS:
!     ==================================================================
      LENGTH = GET_LEN(SIMUNAME)
      
!     PATHS:      
      OUT_SOA     = 'outputs/out.SOA.'//SIMUNAME(1:LENGTH)
      OUT_NCONC   = 'outputs/out.NCONC.'//SIMUNAME(1:LENGTH)
      OUT_PREC    = 'outputs/out.PREC.'//SIMUNAME(1:LENGTH)
      OUT_MDIST   = 'outputs/out.MDIST.'//SIMUNAME(1:LENGTH)
      OUT_NDIST   = 'outputs/out.NDIST.'//SIMUNAME(1:LENGTH)
      OUT_SOMGC   = 'outputs/out.SOMGC.'//SIMUNAME(1:LENGTH)
      OUT_DB      = 'outputs/out.DB.'//SIMUNAME(1:LENGTH)
      OUT_AEMASS  = 'outputs/out.AEMASS.'//SIMUNAME(1:LENGTH)
      OUT_GASDIME = 'outputs/out.GASDIME.'//SIMUNAME(1:LENGTH)
      OUT_GCMASS  = 'outputs/out.GCMASS.'//SIMUNAME(1:LENGTH)
      
!     INITIALIZE OUTFILIES:
!     ==================================================================      
      OPEN(UNIT=40,FILE=OUT_SOA)
      WRITE(40,'(7(A,4X))') 'time','SOA','SOA-W','O2C','ONO','CNO','OLIG'
      CLOSE(40)
      
      OPEN(UNIT=40,FILE=OUT_NCONC)
      WRITE(40,'(3(A,4X))') 'time','NUM_SUSP','NUM_WALL'
      CLOSE(40)

      OPEN(UNIT=40,FILE=OUT_PREC)
      WRITE(40,'(2(A,4X))') 'time','PREC'
      CLOSE(40)

      OPEN(UNIT=44,FILE=OUT_MDIST)
      CLOSE(44)

      OPEN(UNIT=44,FILE=OUT_NDIST)
      CLOSE(44)
      
      FMT = '(0000(A,4X))'; WRITE(FMT(2:5),'(I4)') nSOM + 1

      OPEN(UNIT=40,FILE=OUT_SOMGC)
      WRITE(40,FMT) 'time', (SOMNAMES(j),j=1,nSOM)
      CLOSE(40)
      
      OPEN(UNIT=40,FILE=OUT_DB)
      WRITE(40,'(2(A,4X))') 'time','DB'
      CLOSE(40)
      
      WRITE(FMT(2:5),'(I4)') iCOMP + 1

      OPEN(UNIT=40,FILE=OUT_AEMASS)
      WRITE(40,FMT) 'time', (TOMASNAMES(j),j=1,iCOMP)
      CLOSE(40)
      
      WRITE(FMT(2:5),'(I4)') iCOMP

      OPEN(UNIT=40,FILE=OUT_GCMASS)
      WRITE(40,FMT) 'time', (TOMASNAMES(j),j=1,iCOMP-1)
      CLOSE(40)
      
      END SUBROUTINE
