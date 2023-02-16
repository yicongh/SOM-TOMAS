C     =======================================================================
C     This subroutine determines the condensational driving force of organics
C     between gas and aerosol phases. It then calls the algorithm for
C     condensation (SOLVE_COND).
C      
C     WRITTEN BY Jeff Pierce, Sep. 2010
C     =======================================================================
      
      SUBROUTINE STEP_SOACOND(Pres,Temp,RH,BoxVol,Alpha,Db,dt)

      USE hTOMAS; USE hSOM; IMPLICIT NONE

      DOUBLE PRECISION Pres
      DOUBLE PRECISION Temp
      DOUBLE PRECISION RH
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION Alpha
      DOUBLE PRECISION Db
      DOUBLE PRECISION dt
      
      DOUBLE PRECISION Nki(iBINS), Mki(iBINS,iCOMP), Gci(iCOMP-1)
      DOUBLE PRECISION Nkf(iBINS), Mkf(iBINS,iCOMP), Gcf(iCOMP-1)


      
      DOUBLE PRECISION storg
      PARAMETER(storg=0.025)
      

      double precision dp(ibins)  !Driving force for condensation (Pa)
      double precision tau(ibins)          !condensation parameter (see cond.f)
      double precision atau(ibins)  !same as tau, but all species
      double precision atauc(ibins) !same as atau, but for const dp
      double precision time                !amount of time (s) that has been simulated
      double precision cdt                 !internal, adaptive time step
      double precision mu                  !viscosity of air (kg/m s)
      double precision mfp                 !mean free path of air molecule (m)
      double precision Kn                  !Knudsen number of particle
      double precision Dpk(ibins)          !diameter (m) of particles in bin k
      double precision Rpk(ibins)          !Radius (m) of particles in bin k
      double precision density             !density (kg/m3) of particles
      integer i,k,j,jj,jjj,jo,zz           !counters
      double precision tj, tk(ibins)  !factors used for calculating tau
      double precision sK      !exponential decay const for gas species(g)
      double precision pi, R   !constants
      double precision zeta13  !from Eqn B30 of Tzivion et al. (1989)
      double precision Di, DGAS !diffusivity of gas in air (m2/s)
Ckpc  expand density to all specie
!      double precision gmw(icomp-idiag)        !molecular weight of condensing gas
!      double precision Sv(icomp-idiag)         !parameter used for estimating diffusivity
      double precision beta(ibins)              !correction for non-continuum
      double precision mp,mps(ibins) !particle mass (kg)
      double precision Nko(ibins), Mko(ibins, icomp), Gco(icomp-idiag) !output of cond routine
      double precision mi, mf  !initial and final aerosol masses (updates Gc)
      double precision tr      ! used to calculate time step limits
      double precision mc, ttr
      double precision Neps     !value below which Nk is insignificant
      double precision cthresh  !determines minimum gas conc. for cond.
cdbg      character*12 limit        !description of what limits time step
      double precision tdt      !the value 2/3
      double precision Ntotf, Ntoto, dNerr  !used to track number cons.
      double precision Mktot !Ckpc: total mass for all species           
cdbg      integer numcalls          !number of times cond routine is called
      double precision zeros(ibins) ! array of zeros
      double precision tot_n
      double precision CS       !condensation sink [s-1]
      double precision betat(ibins)
      double precision Mkdry    !dry mass
      double precision WR(ibins) !ratio of wet mass to dry mass of particle
      double precision orgmass ! total mass of organics
      double precision totphase ! total amount of mass for organics to condense into
      double precision pamb ! ambient vapor pressure of species (Pa)
      double precision psat ! saturation vapor pressure of species (Pa)
      double precision totmass, Gcfs
      double precision madd(ibins) ! mass to add or remove from particle
      double precision scalefactor ! scale factor for Psat change by surface tension
      double precision sinkfrac(ibins) ! fraction of condensation sink from each bin
      double precision ms ! mean speed [m/s]
      double precision masseqm ! mass of species in particle at eqm
      double precision masscheqm ! change in mass of species to
                                        ! reach eqm
      double precision maddEQ(ibins) ! change in mass of species per particle
                              ! to reach eqm
      double precision fodc ! first order decay constant towards eqm
	  
      double precision tau_p ! report tau, added by Emily
      double precision test_time
      double precision kB, Na, mair, dorg, dair, morg(iMONO) !(iorg)

      PARAMETER(kB = 1.38E-23)
      PARAMETER(Na = 6.023E23) 

      DOUBLE PRECISION ACCOM       ! accomodation coefficient
      DOUBLE PRECISION FC(IBINS)   ! Fuch's correction factor
      DOUBLE PRECISION KGAS(IBINS)  ! gas-side mass transfer coefficient, m/s
      DOUBLE PRECISION KTOT(IBINS) ! gas-side mass transfer coefficient, m/s
      DOUBLE PRECISION KPAR(IBINS)  !

      
      DOUBLE PRECISION kc          ! first-order loss rate of species in the particle phase, 1/s
      DOUBLE PRECISION qk(IBINS)   ! diffuso-reactive term, unitless
      DOUBLE PRECISION coth(IBINS) ! hyperbolic cotangent
      DOUBLE PRECISION Qkk         ! ratio of particle concentration to surface 
                                   ! concentration at steady state [unitless]

      DOUBLE PRECISION PSAT_ORG(iMONO) !(iORG)
      DOUBLE PRECISION T1
      DOUBLE PRECISION boxmass
      
      DOUBLE PRECISION NonOrgScale
      PARAMETER(NonOrgScale = 0.0)
      
C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------
      double precision CALC_AEDENS
      external CALC_AEDENS

      DOUBLE PRECISION MkORG,MkINORG
      DOUBLE PRECISION m_ORG,m_INORG
      DOUBLE PRECISION Vol_ORG,Vol_INORG,Vol_TOT
      DOUBLE PRECISION r1,r2,dr
      
      DOUBLE PRECISION Vol_COND
      DOUBLE PRECISION Num_COND

      INTEGER CORE_SHELL
      PARAMETER(CORE_SHELL = 1)

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(Neps=1.0d-5, zeta13=0.98483, cthresh=1.d-16)
      
      DOUBLE PRECISION FACTOR1,FACTOR2

      DOUBLE PRECISION Nepss
      PARAMETER(Nepss = 1.e-3)
      
      DOUBLE PRECISION LRATIO
      DOUBLE PRECISION BETAS
      DOUBLE PRECISION ppp
      
      DOUBLE PRECISION,EXTERNAL :: CALC_KPAR
      DOUBLE PRECISION,EXTERNAL :: CALC_QBULK

      DOUBLE PRECISION u1,u2
      DOUBLE PRECISION mtotal

      DOUBLE PRECISION kcc
C      PARAMETER(kcc = 0.01)

      DOUBLE PRECISION QQQ(iBINS)

      DOUBLE PRECISION Csat_Org(iMONO) !(iORG)
      DOUBLE PRECISION Mass_Org
      
      DOUBLE PRECISION kkkc
      DOUBLE PRECISION xxxx

      DOUBLE PRECISION QBULK(iBINS)
      DOUBLE PRECISION QFACTOR

      DOUBLE PRECISION OUT1,OUT2
      DOUBLE PRECISION SRATIO

      INTEGER xC10O4

      INTEGER MASK(iBINS)

      DOUBLE PRECISION Rp,Rc

      INTEGER INDEX
      
      
      do k=1,ibins
         if (Nk(k) .lt. Nepss) then
            Nk(k)=Nepss
            do j=1,icomp
               Mk(k,j)=0.d0
            enddo
            Mk(k,xSO4)=Nepss*1.4*xk(k) !make the added particles SO4 !xSO4 is an index that points to the position of SO4
         endif                          !xSO4 = 2
      enddo

      do k=1,ibins
         do j=1,icomp-idiag
            if (Mk(k,j) .eq. 0.d0) then
               !add a small amount of mass
               mtotal=0.d0
               do jj=1,icomp-idiag
                  mtotal=mtotal+Mk(k,jj)
               enddo
               Mk(k,j)=1.d-10*mtotal
            endif
         enddo
      enddo
      
C-----CODE--------------------------------------------------------------

      Mass_ORG = SUM(Mk(:,xORG_1ST:xORG_LST))
      

      T1 = 298.0

      kcc = kc
      
      boxmass = 0.0289 * pres * BOXVOL * 1.0e-6 / R /temp


      
 201  CONTINUE

      
      DO j = 1,iMONO
         PSAT_ORG(j) = PSAT_REF(j)*EXP(Hvap(j)*1.d3*(1.0/T1 - 1.0/Temp)/R)
         CSAT_ORG(j) = CSAT_REF(j)
      END DO

      
      

      dNerr=0.0
      tdt=2.d0/3.d0

C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gcf(j)=Gc(j)
      enddo


      
      do k=1,ibins
         Nkf(k)=Nk(k)
         do j=1,icomp
            Mkf(k,j)=Mk(k,j)
         enddo
      enddo

      

C Calculate tj and tk factors needed to calculate tau values
C      mu=2.5277e-7*temp**0.75302 !original
       mu = 1.458e-6*temp**1.5/(110.4 + temp)
      !mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
       mair = (0.79*28. + 0.21*32.)/1000./Na  ! molecular weight of air, kg/mol, Andy
        
       dorg = 10e-10            ! diameter of organic molecule, m,Andy
       dair = (0.79*1.09e-10 + 0.21*1.21e-10) ! diameter of air molecule, m,Andy	  

      
      do jo=1,iMONO !iorg-2
      
      time=0.0  !go to next species when time=dt
C      ORG_MW(jo)=(0.434 - 0.045*log10(cstar(jo)))*1000. !keep consistent with wall loss module
      morg(jo) = ORG_MW(jo)/Na/1000.0		  !keep consistent with wall loss module

      j=xORG_1ST+jo-1
      ms=sqrt(8.0*R*temp/(pi*ORG_MW(jo)*1.0d-3)) !S&P2 eqn 9.2
C      Di=diorg(jo)*(temp/T1)**ditemporg(jo)      !J's diffusivity of organic vapors
C	  Di = (2./3.)*sqrt(kB*temp/pi*0.5*(1./morg(jo)+1/mair))*1./ 
C     &	       (pi*(0.5*(dorg + dair))**2.)/Na*(R*temp/pres) ! Andy's diffusivity of organic vapors, m2/s
      Di = 1.38e-5*44.0/ORG_MW(jo)

      
      mfp=3.d0*Di/ms
      
      tj=2.*pi*ORG_MW(jo)*1.0d-3/(R*temp)
      
      
C Repeat from this point if multiple internal time steps are needed
 10   continue

      CALL MNFIX(Nkf,Mkf)
      
      sK=0.0d0
      


C     OVERALL MASS TRANSFER COEFFICIENT:
C     =======================================================================
      DO k = 1,iBINS
	     
         if (Nkf(k) .gt. Neps) then

Ckpc  calculate the total mass   
Ckpc   factor of 1.2 assumes ammonium bisulfate.
Ckpc  Add 0.2x first, and then add 1.0x in the loop
            Mktot=0.d0
            Mkdry=0.d0
            do jj=1,icomp
               Mktot=Mktot+Mkf(k,jj)
            enddo
            
            do jj=1,icomp-idiag
               Mkdry=Mkdry+Mkf(k,jj)
            enddo
            WR(k)=Mktot/Mkdry
             
            orgmass=0.d0
            
            do jj=1,iorg
               orgmass=orgmass+Mkf(k,xORG_1ST+jj-1) 
            enddo
            
C            density=aerodens(Mkf(k,xSO4)+orgmass,0.d0,
C     &         Mkf(k,xNH4),0.d0,  !assume bisulfate and org is sulf acid
C     &         Mkf(k,xH2O))

            density = CALC_AEDENS(orgmass,Mk(k,1))
 
            mp=Mktot/Nkf(k)

            
            
         else
            !nothing in this bin - set to "typical value"

C     THIS VALUE HAS BEEN CHANGED!
C            density=1500.
            density=1700.
            
C            mp=1.4*xk(k)
C            mp=1.65*xk(k)
            mp = SQRT(xk(k)*xk(k+1))
            
         endif
         
         mps(k)=mp
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
C         print*, '-------------------------'
C         print*, 'Dpk(k) =', Dpk(k)
C         print*, '-------------------------'
         Kn=2.0*mfp/Dpk(k)      !S&P eqn 11.35 (text)

         

!     AliA - calculate KGAS
CCCCC ============================================================================================
!     basic properties; can be eliminated if defined in TOMAS prior to this
         Rpk(k) = Dpk(k)/2.0    ! radius, m
         DGAS = Di               ! gas-phase diffusion coefficient, m2/s (constant? if yes, can we vary this based on information about the SOM species)
         accom = alpha          ! mass accommodation coefficient, unitless


         FC(k) = 0.75 * accom * (1.0+Kn) / (Kn*(1.0+Kn)
     &        + 0.238*accom*Kn + 0.75*accom) ! Fuchs correction, unitless, eqn (14)

         
C        CALCULATE PARTICLE AND CORE RADII:
         Vol_ORG   = SUM(Mk(k,xORG_1ST:xORG_LST))/Nk(k)/1400.
         Vol_INORG = Mk(k,1)/Nk(k)/1770.
         Vol_TOT   = Vol_ORG + Vol_INORG
         
         Rp = (3.0*Vol_TOT/4.0/PI)**(1./3.)
         Rc = (3.0*Vol_INORG/4.0/PI)**(1./3.)



C        FIRST-ORDER DECAY RATE:         
         kc = kc_HET(k,xORG_1ST+jo-1) + kc_OLIG(k,xORG_1ST+jo-1)
         
C        GAS-SIDE MASS-TRANSFER COEFFICIENT:         
         KGAS(k) = DGAS/Rpk(k)*FC(k)
         
C        PARTICLE-SIDE MASS-TRANSFER COEFFICIENT: 
         KPAR(k)  = CALC_KPAR(Rc,Rp,Db,kc)
         QBULK(k) = CALC_QBULK(Rc,Rp,Db,kc)
         
C        OVERALL MASS-TRANSFER COEFFICIENT:
         SRATIO  = Csat_Org(jo)/density/1.0e9*QBULK(k)
         KTOT(k) = Rp/(1.0/KGAS(k) + 1.0/KPAR(k)*SRATIO)
         
                  
CCCCC ============================================================================================
!     AliA - calculate KGAS
         
Ckpc  adjust beta, tk, sK and alpha for all species
C         beta(k)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha) !S&P eqn 11.35
         
         tk(k)=(6./(pi*density))**(1./3.)*KTOT(k)
         
         if (Nkf(k) .gt. 0.0) then
            sK=sK+Dpk(k)*Nkf(k)*KTOT(k)
         endif
         
      enddo
      
      do k=1,ibins
!         sinkfrac(k)=tk(k)*Nkf(k)*mps(k)**(1.d0/3.d0)/sK
         sinkfrac(k)=Dpk(k)*Nkf(k)*KTOT(k)/sK
      enddo
      
!      sK=sK*zeta13*tj*R*temp/(ORG_MW(jo)*1.d-3)/(boxvol*1.d-6)
      sK = 2.d0*pi*sK/(dble(boxvol)*1D-6)

      

C Set dp for nonvolatile species
      DO k = 1,iBINS

         TOTPHASE = SUM(Mk(k,xORG_1ST:xORG_LST)) + 
     &              Mk(k,xSO4)*NONORGSCALE
         
         
         scalefactor=exp((storg*0.001*ORG_MW(j))
     &                     /(8.314*temp*1200.*Dpk(k))) !! kelvin effect, mw dependent
         
         pamb=(Gcf(j)/ORG_MW(jo))/(boxmass/28.9)*pres 
         
         
         if (totphase.gt.0.d0) then
            psat=PSAT_ORG(jo)*Mkf(k,j)/totphase*scalefactor/QBULK(k) ! This is Raoult's law
C            psat=PSAT_ORG(jo)*(Mkf(k,j)+Mkf(k,xDIME_1ST+jo-1))/
C     &      totphase*scalefactor ! This is Raoult's law
C            PRINT*, 'MONO=',Mkf(k,j),'DIMER=',Mk(k,xDIME_1ST+jo-1)
C            psat=PSAT_ORG(jo)*(Mkf(k,j)+Mk(k,xDIME_1ST+jo-1))
C     &                /totphase*scalefactor ! This is Raoult's law

            !psat=PSAT_ORG(iORG-2)*Mkf(k,j)/totphase*scalefactor ! This is Raoult's law
            
         else
            psat=0.d0
         endif
         
         dp(k) = pamb - psat !!!!!! NOTE CHANGE BY CHARLES
         
         !!! do not let evaporate for now !!!
         !!! There are big issues with evaporation for multiple species with multiple size bins!
         !!! I need to find a way to limit the timestep so that species don't entirely evaporate,
         !!! and the system does not get overly stiff
!         if (dp(k).lt.0.d0) then
!            dp(k) = 0.d0
!         endif
         
         ! Determine if equillibrium amount of mass of this species exists 
         ! (holding all else constant).  What mass is it?
         ! For now assume surface tension scalefactor is constant
         ! Equillibrium occurs when dp(k) = 0
        if (totphase .EQ. 0.0) then
		 masseqm = 0.0
		 masscheqm = 0.0
		 maddEQ(k) = 0.0
       else 

         
         masseqm=pamb*QBULK(k)/(PSAT_ORG(jo)/totphase*scalefactor)
C         masseqm=pamb/(PSAT_ORG(jo)/totphase*scalefactor)
          
         !masseqm=pamb/(PSAT_ORG(iORG-2)/totphase*scalefactor)
         
         masscheqm=masseqm-Mkf(k,j) ! change in mass to reach eqm
         maddEQ(k)=masscheqm/Nkf(k)  !!!!!CHANGED BY CHARLES
       endif
!         print*, j,k,dp(k),Mkf(k,j),masseqm,masscheqm,maddEQ(k),Nkf(k)
!         if ((jo.eq.5).and.(k.eq.15)) then
!		    print*,'test',k,pamb,psat,masseqm,Mkf(k,j)/Nkf(k),maddEQ(k)
!     endif
       
      enddo
C      STOP

      


      cdt = dt
C      cdt = dt/2.0
C      cdt = dt/10.0

      do k=1,ibins
         
         !!! NAN CHECK !!!
         if (isnan(sK)) then
            print*,'TAU sK in soacond'
            print*,'j=',j
            print*,'time=',time
            PRINT*,'Kg=',KTOT(k),'kg=',KGAS(k),'kp=',KPAR(k),'Q=',QQQ(k)
            print*,'kc=',kc,'sK=',sK,'NumCond=',Num_COND
            PRINT*,'Vol_COND=',Vol_COND,'Nk=',Nk(k)
            PRINT*, 'Vol_ORG=',Vol_ORG,'Mk=',Mk(k,:)
C            print*,'tk=',tk
C            print*,'tj=',tj
C            print*,'dp=',dp
C            print*,'cdt=',cdt
C            print*,'WR=',WR
C            print*,'Nk=',Nk
C            print*,'Mk=',Mk
            stop
         endif

         
         if (sK .gt. 0.0) then
               atau(k)=tj*tk(k)*dp(k)/sK
     &             *(1.d0-exp(-1.d0*sK*cdt))*WR(k)
         else
           !!nothing to condense onto
            atau(k)=0.0
         endif


         
         mc=0.0d0
         do jj=1,icomp
            mc=mc+Mkf(k,jj)/Nkf(k)
         enddo
         
         
         
         if (atau(k)/1.5d0+mc**tdt .ge. 0.d0) then
            madd(k)=((atau(k)/1.5d0+mc**tdt)**(1.d0/tdt)-mc)/WR(k)
         else
C            madd(k)=((atau(k)/1.5d0+mc**tdt)**(1.d0/tdt)-mc)/WR(k)            
            madd(k)=maddEQ(k) !!!!!! EDITED BY CHARLES
         endif
         

         if (Nkf(k).lt.Neps)then
            atau(k)=0.d0
            madd(k)=0.d0
         else
            if((((maddEQ(k).gt.0.d0).and.(madd(k)).lt.0.d0).or.
     &       ((maddEQ(k).lt.0.d0).and.(madd(k).gt.0.d0))).and.
     &       (dp(k)/pamb.gt.1.d-4))then
C               print*,'maddEQ and madd have opposite signs in 
C     &h2oconda'
C               print*,'k',k,'j',j
C               print*,'maddEQ',maddEQ(k),'madd',madd(k)
C               print*,'dp',dp(k),'pamb',pamb
c               STOP
               xxxx=1.
            endif

C            PRINT*, k,maddEQ(k),madd(k)
            
            if (madd(k)/maddEQ(k).gt.1E-30)then
                                                     ! 1st order model
                                                     ! 2nd condition
                                                     ! tests for NAN
               fodc=1.d0/(cdt*(maddEQ(k)/madd(k))) ! first
                   ! order decay constant towards equil
               madd(k)=maddEQ(k)*(1.d0-exp(-fodc*cdt)) ! change
                         ! to new mass to add
c               print*,'k',k,'atau1',atau(k,xH2O)
               atau(k)=1.5d0*((mc+madd(k)*WR(k))
     &              **(2.d0/3.d0)-mc**(2.d0/3.d0)) ! find new tau
                             !value
c               print*,'k',k,'atau2',atau(k,xH2O)
            endif
         endif
         
!         print*,'k',k,'madd',madd(k),'atau',atau(k)

!         if (-madd(k).gt.Mkf(k,j)/Nkf(k)*1.000001)then  !original 
          if ((-madd(k)-Mkf(k,j)/Nkf(k)) .GT. 1.00d-17) then !modified by Emily
            print*,'too much evaporation in SOACOND!'
            print*,'k',k,'j',j
            print*,'madd',madd(k)
            print*,'mp',Mkf(k,j)/Nkf(k)
            print*,'Mk',Mkf(k,j),'NKf',Nkf(k)
            PRINT*, QBULK(k)
            stop
         endif

C         PRINT*, k,atau(k)
         
      enddo
C      STOP

C Call condensation subroutine to do mass transfer


       !Swap tau values for this species into array for cond
      do k=1,ibins
         tau(k)=atau(k)
         !!! NAN CHECK !!!
         if (isnan(tau(k))) then
            print*,'TAU NaN in soacond'
            print*,'j=',j,'k=',k
            print*,'time=',time
            print*,'tau=',atau
            print*,'sK=',sK
            print*,'tk=',tk
            print*,'tj=',tj
            print*,'dp=',dp
            print*,'cdt=',cdt
            print*,'WR=',WR
            print*,'Nk=',Nk
            print*,'Mk=',Mk
            stop
         endif
      enddo


      !Call condensation routine
      Ntotf=0.0
      do k=1,ibins
         Ntotf=Ntotf+Nkf(k)
         !!! NAN CHECK !!!
         if (isnan(Nkf(k))) then
            print*,'Nk NaN in soacond before tmcond'
            print*,'j=',j,'k=',k
            print*,'time=',time
            print*,'tau=',atau
            print*,'sK=',sK
            print*,'tk=',tk
            print*,'tj=',tj
            print*,'dp=',dp
            print*,'cdt=',cdt
            print*,'WR=',WR
            print*,'Nkf=',Nkf
            print*,'Mkf=',Mkf
            stop
         endif
      enddo


C     CONDENSATION MASK:
      MASK(:) = 1


      
      CALL TMCOND(TAU,Xk,Mkf,Nkf,Mko,Nko,j,MASK)


      !Check for number conservation
      Ntoto=0.0
      do k=1,ibins
         Ntoto=Ntoto+Nko(k)
         !!! NAN CHECK !!!
         if (isnan(Nkf(k))) then
            print*,'Nk NaN in soacond after tmcond'
            print*,'j=',j,'k=',k
            print*,'time=',time
            print*,'tau=',atau
            print*,'sK=',sK
            print*,'tk=',tk
            print*,'tj=',tj
            print*,'dp=',dp
            print*,'cdt=',cdt
            print*,'WR=',WR
            print*,'Nkf=',Nkf
            print*,'Mkf=',Mkf
            stop
         endif
      enddo
         
      dNerr=dNerr+Ntotf-Ntoto
      
      if (abs(dNerr/Ntoto) .gt. 1.e-4) then
         write(*,*) 'ERROR in soacond: Number not conserved'
         write(*,*) 'time=',time
         write(*,*) 'j=', j
         write(*,*) Ntoto, Ntotf
         write(*,*) (Nkf(k),k=1,ibins)
      endif

      do k=1,ibins
         if (Nko(k).lt.0.d0) then
            print*, 'Nk < 0 in soacond'
            print*, 'k=',k,'j=',j
            print*,'Nk=',Nko(k)
            print*,'Mk=',Mko(k,:)
            stop
         endif
         do jj=1,icomp
            if (Mko(k,jj).lt.-1.0d-10) then
               print*, 'Mk < 0 in soacond'
               print*, 'k=',k,'j=',j,'jj=',jj
               print*,'Nkf=',Nkf(k)
               print*,'Mkf=',Mkf(k,:)
               print*,'Nko=',Nko(k)
               print*,'Mko=',Mko(k,:)
               mc=0.0d0
               do jjj=1,icomp
                  mc=mc+Mkf(k,jjj)/Nkf(k)
               enddo
               print*,'mc',mc
               print*,'mp',Mkf(k,jj)/Nkf(k)
               print*,'madd1',madd(k)
               print*,'madd2',((atau(k)/1.5d0+mc**tdt)**(1.d0/tdt)-mc)
     &                         /WR(k)
               print*,'time=',time
               print*,'tau=',atau(k)
               print*,'sK=',sK
               print*,'tk=',tk
               print*,'tj=',tj
               print*,'dp=',dp
               print*,'cdt=',cdt
               print*,'WR=',WR
               stop
               !Mko(k,jj)=0.d0
            endif
         enddo
      enddo

      !Update gas phase concentration
      mi=0.0
      mf=0.0
      do k=1,ibins
         mi=mi+Mkf(k,j)
         mf=mf+Mko(k,j)
      enddo
      Gcfs=Gcf(j)
      Gcf(j)=Gcf(j)+(mi-mf)
      if (Gcf(j).lt.0.d0)then
         if ((abs(Gcf(j)/mi).lt.1d-8).or.(abs(Gcf(j)).lt.1d-8))then
            Gcf(j)=0.d0
         else
            print*,'soacond Gc < 0'
            print*,'j',j
            print*,Gcfs,Gcf(j),(mi-mf),mi,mf
            print*,'cdt',cdt
            print*,'tau',tau
            stop
         endif
      endif
      
      
      !Swap into Nkf, Mkf
      do k=1,ibins
         Nkf(k)=Nko(k)
         do jj=1,icomp
            Mkf(k,jj)=Mko(k,jj)
         enddo
      enddo
      

      
      
      !Update ammonia concentrations
C      call eznh3eqm(Gcf,Mkf)

      !Update water concentrations
C      call ezwatereqm(Mkf)

C Update time
      time=time+cdt
      
C Repeat process if necessary
C      if (time .lt. dt) goto 10

      test_time = time - dt
      if (abs(test_time).lt.1.00d-8) then
        test_time = 0.00d0
      endif
 
       if (test_time.lt.0.00d0) goto 10

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo ! loop over components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
cdbg      write(*,*) 'Cond routine called ',numcalls,' times'
cdbg      write(*,*) 'Number cons. error was ', dNerr

 100  continue   !skip to here if there is no gas phase to condense

      do k=1,ibins
         Nk(k)=Nkf(k)
         do jj=1,icomp
            Mk(k,jj)=Mkf(k,jj)
         enddo
      enddo
      
      Gc(:) = Gcf(:)
      
      RETURN
      END SUBROUTINE






