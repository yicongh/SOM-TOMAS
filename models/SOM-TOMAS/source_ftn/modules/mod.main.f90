!     ==========================================================================
!                               THIS IS THE MAIN MODULE
!     ==========================================================================

      MODULE hMAIN
      
      USE hINVARS
      USE hSOM
      USE hTOMAS
      
      IMPLICIT NONE

      REAL(8) :: odt
      REAL(8) :: timer

      REAL(8) :: OH
      REAL(8) :: JNUC

      INTEGER :: iCTR
      INTEGER :: i,j,k

      REAL(8),EXTERNAL :: CALC_JNUC

      REAL(8),PARAMETER :: NA = 6.022e23
      REAL(8),PARAMETER :: R = 8.314

      END MODULE hMAIN
