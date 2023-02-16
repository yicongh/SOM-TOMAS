!     =======================================================================
!                  THIS SUBROUTINE STEPS FORWARD FOR NUCLEATION
!     =======================================================================
      
      SUBROUTINE STEP_NUCL(JNUC,BoxVol,timer,mdt)

      USE hTOMAS; IMPLICIT NONE
      
      DOUBLE PRECISION JNUC
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION timer
      DOUBLE PRECISION mdt

      DOUBLE PRECISION dN
      
      dN = JNUC*BoxVol*mdt
      
      Nk(1) = &
      Nk(1) + dN

      Mk(1,xORG_LST) = & 
      Mk(1,xORG_LST) + dN*SQRT(Xk(1)*Xk(2)) ! NOTE THE INDEX HERE - CHECK ITS C*
      
      RETURN
      END SUBROUTINE
