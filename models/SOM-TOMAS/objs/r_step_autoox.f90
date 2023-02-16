!     =======================================================================
!       THIS SUBROUTINE STEPS FORWARD FOR AUTO-OXIDATION OF THE PRECURSORS
!     =======================================================================
      
      SUBROUTINE STEP_AUTOOX(f_HOMS,PRECS_kOH,OH,Temp,Pres,BoxVol,timer,dt)

      USE hTOMAS; USE hSOM; IMPLICIT NONE
      
      DOUBLE PRECISION f_HOMS(nSOMPRECS)
      DOUBLE PRECISION PRECS_kOH(nSOMPRECS)
      DOUBLE PRECISION dt
      DOUBLE PRECISION OH
      DOUBLE PRECISION Temp
      DOUBLE PRECISION Pres
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION timer
            
      DOUBLE PRECISION NA
      PARAMETER(NA = 6.02E23)

      DOUBLE PRECISION R
      PARAMETER(R = 8.314)
      
      INTEGER i,j,k

      DOUBLE PRECISION PREC,dPREC
      INTEGER INDEX,xHOMS(nSOMCLASS)
      INTEGER INDEX1,INDEX2

      DOUBLE PRECISION dSO2

      IF (timer.LT.0.0) GOTO 100
      
!     STEP FOR AUTO-OXIDATION:
!     =======================================================================
      k = 1
      
      DO i = 1,nSOMCLASS
         DO j = 1,uSOMPRECS(i)
!           PRECURSOR INDEX:
            INDEX = xPREC_1ST(i)+j-1

!           PRECURSOR REACTED:
            PREC  = SOMGC(INDEX)
            dPREC = PRECS_kOH(k)*7.34e15/300.*PREC*OH*dt*f_HOMS(k) !!!!!!!!

!           DEDUCT FROM PRECURSOR:
            SOMGC(INDEX) = SOMGC(INDEX) - dPREC      
            k = k + 1
         END DO

!        ADD TO HOMS:         
         INDEX1 = xGRID_LST(i)-xGRID_1ST(1)+xORG_1ST
         INDEX2 = xGRID_LST(i)-xGRID_1ST(1)+1

         GC(INDEX1) = GC(INDEX1) + &
         dPREC*1e-6*(Pres/Temp/R)*ORG_MW(INDEX2)*0.001*(BoxVol*1e-6)

      END DO


!     PROXY FOR H2SO4:
      dSO2 = 9.0E-13*7.34e15/300.*SO2*OH*dt  !!!!!!!!!!!!!!!!!
      
      SO2 = &
      SO2 - dSO2
      
      GC(2+iMONO-1-1) = &
      GC(2+iMONO-1-1) + dSO2*1e-6*(Pres/Temp/R)*98.*0.001*(BoxVol*1e-6)

      
 100  CONTINUE
      RETURN
      END SUBROUTINE
