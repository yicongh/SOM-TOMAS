C     ================================================================
C     This subroutine examines the mass and number distributions and
C     determines if any bins have an average mass outside their normal
C     range.

C     WRITTEN BY Peter Adams, September 2000
C     ================================================================
      
      SUBROUTINE MNFIX(Nkx,Mkx)

      USE hTOMAS; IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      DOUBLE PRECISION Nkx(iBINS)
      DOUBLE PRECISION Mkx(iBINS,iCOMP)
      DOUBLE PRECISION Okx(iORG,iORG,iBINS)


      integer i,j,k,n,kk           !counters
      integer newbin           !bin number into which mass is shifted
      double precision Xold, Xnew  !average masses of old and new bins
      double precision DryMass !dry mass of in a bin
      double precision avg      !average dry mass of particles in bin
      DOUBLE PRECISION AvgMass
      double precision number  !number of particles initially in problem bin
      double precision nshift  !number to shift to new bin
      double precision mshift  !mass to shift to new bin
      double precision fj      !fraction of mass that is component j
      double precision eps     !small number
      double precision Neps

      DOUBLE PRECISION PRATIO
      DOUBLE PRECISION MassFrac
      
      parameter(eps=1.d-40,Neps=1E-5)

C === BIN-MULTIPLYING RATIO ==========================================

      PRATIO = Xk(2)/Xk(1)

CCCCC NOTE HERE !!!!!!!!!!!!!!!!!!!!!!!1
C      RETURN
      
      
      
C === IF NUMBER IS TINY IN BINS ======================================
      
      DO i = 1,iBINS
         IF (Nkx(i).LT.Neps) THEN
            
            DO j = 1,xORG_LST
               IF (Mkx(i,j).LT.1E-5) THEN                  
                  Mkx(i,j) = 0.0
               ELSE
                  WRITE(*,*) 'N=0 but M>1.d-5 in mnfix'
                  WRITE(*,*) 'bin=',i
                  STOP
               END IF
            END DO
            
            Nkx(i)      = Neps
            Mkx(i,xSO4) = Neps*SQRT(Xk(i)*Xk(i+1))
            
         END IF
      END DO

C === IF AVE. BIN MASS IS COMPLETELY OUT OF RANGE ====================
      
      DO i = 1,iBINS

         DryMass = SUM(Mkx(i,1:xORG_LST)) +
     &             SUM(Okx(:,:,i))    
         AvgMass = DryMass/(Nkx(i) + eps)
         
         ! IF OUT OF THE UPPER BIN RANGE (REMOVE SOME MASS):
         IF (AvgMass.GT.Xk(iBINS+1)) THEN
            
            MassFrac = Nkx(i)*Xk(iBINS+1)/(0.9 + 0.1*PRATIO)/
     &                 (DryMass + eps)
            
            DO j = 1,iCOMP
               Mkx(i,j) = Mkx(i,j) * MassFrac
            END DO
            
         END IF

         ! IF OUT OF THE LOWER BIN RANGE (REMOVE SOME NUMBER):
         IF (AvgMass.LT.Xk(1)) THEN
            
            Nkx(i) = DryMass/(Xk(1)*(0.9 + 0.1*PRATIO))
            
         END IF
         
      END DO

C === IF AVE. BIN MASS IS OUT OF BIN BOUNDS ==========================
      
      DO i = 1,iBINS

         DryMass = SUM(Mkx(i,1:xORG_LST)) !+
C     &             SUM(Okx(:,:,i))    
         AvgMass = DryMass/(Nkx(i) + eps)
         
         IF (Nkx(i).EQ.0.0) THEN
            IF (DryMass.GT.0.0) THEN
               WRITE(*,*) 'N=0 but M>0 in mnfix'
               STOP
            ELSE
               AvgMass = SQRT(Xk(i)*Xk(i+1))
            END IF
         END IF

         ! IF AVERAE MASS IS GREATER THAN UPPER BIN BOUND:
         IF (AvgMass.GT.Xk(i+1)) THEN
            
            n    = i + 1
            Xnew = Xk(n+1)/(0.9 + 0.1*PRATIO)
            
 100        IF (Xnew.LT.AvgMass) THEN
               n    = n + 1
               Xnew = Xk(n+1)/(0.9 + 0.1*PRATIO)
               GOTO 100
            END IF
            
            Xold = SQRT(Xk(i)*Xk(i+1))
            
            nshift = (DryMass - Xold*Nkx(i))/(Xnew - Xold)
            Nkx(i) = Nkx(i) - nshift
            Nkx(n) = Nkx(n) + nshift
            
            mshift = Xnew*nshift
            
            DO j = 1,iCOMP
               fj       = Mkx(i,j)/(DryMass+eps)
               Mkx(i,j) = Xold*Nkx(i)*fj
               Mkx(n,j) = Mkx(n,j) + mshift*fj
            END DO
            
         END IF

         ! IF AVERAE MASS IS LESS THAN LOWER BIN BOUND:
         IF (AvgMass.LT.Xk(i)) THEN
            
            n    = i - 1
            Xnew = Xk(n)*(0.9 + 0.1*PRATIO)
            
 200        IF (Xnew.GE.AvgMass) THEN
               n    = n - 1
               Xnew = Xk(n)*(0.9 + 0.1*PRATIO)
               GOTO 200
            END IF
            
            Xold = SQRT(Xk(i)*Xk(i+1))
            
            number=Nkx(i)
            nshift=(drymass-xold*number)/(xnew-xold)
            mshift=xnew*nshift
            Nkx(i)=Nkx(i)-nshift
            Nkx(n)=Nkx(n)+nshift
            
            do j=1,icomp
               fj=Mkx(i,j)/(drymass+eps)
               Mkx(i,j)=xold*Nkx(i)*fj
               Mkx(n,j)=Mkx(n,j)+mshift*fj
            enddo
            
         endif

         
         
      enddo

      DO k = 1,iBINS
         IF (Nkx(k).LE.1D-5) THEN
            Nkx(k) = 1D-5
         END IF
         DO j = 1,iCOMP
            IF (Mkx(k,j).LE.1D-5*SQRT(Xk(k)*Xk(k+1))) THEN 
               Mkx(k,j) = 1D-5*SQRT(Xk(k)*Xk(k+1))
            END IF
         END DO
      END DO
      
      RETURN
      END
