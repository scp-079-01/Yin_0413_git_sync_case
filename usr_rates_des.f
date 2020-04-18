      SUBROUTINE USR_RATES_DES(NP, pM, IJK, DES_RATES)

      USE compar
      USE constant
      USE des_thermo          !des_T_s des_C_ps
      USE des_rxns
      USE discretelement      !pmass-- mass of one particle
      USE energy
      USE fldvar				! x_s: Solid phase species mass fractions. Also include ROP_S(macroscopic gas density??),RO_s(gas density)	! x_s(IJK,solid phase index, global species index)
      !2020/4/12 x_s is not right, should use des_x_s(NP,spcies)
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE rxns
      USE run
      USE usr
      USE fun_avg
      USE functions

      USE parallel
      USE sendrecv
      USE toleranc
      

      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP  ! Global index of particle
      INTEGER, INTENT(IN) :: pM  ! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP
      
! Calculated reaction rates. (reacted moles per sec)
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)
      
      
! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Solids Diameter (m)
      DOUBLE PRECISION :: Dp0
! Gas phase concentrations (mol/m^3)
      DOUBLE PRECISION c_O2,  c_H2, c_H2O, c_CO2
      DOUBLE PRECISION :: c_Biomass
! Partial pressures (KPa)
      DOUBLE PRECISION p_O2,  p_H2, p_H2O, p_CO2
      
      DOUBLE PRECISION E_i, A_i, C_i, A_p
      DOUBLE PRECISION r_diff, r_kin
      
 ! Bounded phase temperatures (K)
      DOUBLE PRECISION xTg   ! gas phase
      DOUBLE PRECISION xTs   ! solids phase
      
      integer LB,UB

            
! DOUBLE PRECISION :: TsX   ! solids phase 1

! Prevent reactions when mass fraction falls below c_Limiter.
      DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-7

!     Uncomment the following line when customizing this usr_rates_des.f to include generated file species.inc:
!
     INCLUDE 'species.inc' !
     
     
! Reaction rates:
!`````````````````````````````````````````````````````````````````````//
! Calculate particle diameter.

      IF (PM==1) then !Only the biomass would react, sand(PM==2) will not.
      IF(pmass(NP)>1.0d-10) then !prevent too small particles
      
            
      ! Calculate the bounded temperatures.
      !xTg = max(min(T_g(IJK),   TMAX), 300.0d0) ! low bound 300, high bound TMAX=4000(equals to TMAX defined in Cp)
      !xTs = max(min(DES_T_s(NP), TMAX), 300.0d0)
      
      xTg = max(min(T_g(IJK),   2000.0), 300.0d0) 
      xTs = max(min(DES_T_s(NP), 2000.0), 300.0d0)
      
      
      !!2020/4/13 replace all the T_g(IJK) with xTg, DES_T_s(NP) with xTs
      ! 202/4/13 replace all the 12.0 with MW_s(1,Carbon)
     
      Dp0 = 2.0d0 * DES_RADIUS(NP)
      
      write(*,*) ''
      write(*,*) '****************************'
      write(*,*) 'A usr_rates_des is called'
      write(*,*) 'NP= ', NP
      write(*,*) 'PM ', PM
      write(*,*) 'Dp0= ', Dp0
      write(*,*) 'des_T_s(NP)=', des_T_s(NP)
      write(*,*) 'xTg=', xTg
      write(*,*) 'xTs=', xTs
      write(*,*) 'TMAX=', TMAX
            
! Particle surface area (m^2)
      A_p = Pi* Dp0 * Dp0   ! (m^2)
      write(*,*) 'Ap= ', A_p
! Initialize partial pressures and concentrations.
      c_O2   = ZERO;     p_O2 = ZERO;
      c_H2O  = ZERO;     p_H2O = ZERO;
      c_H2   = ZERO;     p_H2   = ZERO;
      c_CO2   = ZERO;    p_CO2   = ZERO;
      

! Calculate the partial pressure and molar concentration of 
      IF(X_g(IJK,H2O) > c_Limiter) then
            p_H2O = 1e-3*P_g(IJK) * X_g(IJK,H2O) * (MW_MIX_g(IJK) / MW_g(H2O))
            c_H2O = RO_g(IJK) * X_g(IJK,H2O) / Mw_g(H2O)
      ENDIF
      write(*,*) 'p_H2O= ', p_H2O
      write(*,*) 'c_H2O= ', c_H2O
      write(*,*) 'X_g(IJK,H2O)= ', X_g(IJK,H2O)
            
! Calculate the partial pressure and molar concentration of 
      IF(X_g(IJK,CO2) > c_Limiter) then
         p_CO2 = 1e-3*P_g(IJK) * X_g(IJK,CO2) * (MW_MIX_g(IJK) / MW_g(CO2))
         c_CO2 = RO_g(IJK) * X_g(IJK,CO2) / Mw_g(CO2)
      ENDIF
      write(*,*) 'p_CO2= ', p_CO2
      write(*,*) 'c_CO2= ', c_CO2
      write(*,*) 'X_g(IJK,CO2)= ', X_g(IJK,CO2)
      
      !c_Biomass = ROP_s(IJK,1)*X_s(IJK,1,BiomassCore)/MW_s(1,BiomassCore)
!**********************************************************************!
!                        Heterogeneous Reactions                       !
!**********************************************************************!
   if(.not. compare(EP_s(IJK,1),ZERO)) then  !only react when there is biomass particle

      write(*,*) 'PM==1'
      write(*,*) 'MW_s(1,Carbon)=', MW_s(1,Carbon)
      write(*,*) 'MW_s(1,2)=', MW_s(1,2)
      write(*,*) 'MW_s(1,3)=', MW_s(1,3)! 1 represents solid phase, 3 represents species:biomassCore
      write(*,*) 'EP_s(IJK,1)=', EP_s(IJK,1)
      write(*,*) 'BiomassCore=', BiomassCore
      write(*,*) 'IJK', IJK
      
! C + CO2 --> CO  (kg-mole/m^3.s)
!---------------------------------------------------------------------//
      if (DES_X_s(NP,Carbon) > c_Limiter) then  !2020/4/17 added  

         C_i=5d-12
         A_i=8.3
         E_i=4.37d7
         r_kin   = A_i *xTs * exp(-E_i/(8.3147295*1000*xTs))
         r_diff  = C_i * ((xTs+xTg)/2)**0.75/Dp0
         DES_RATES(Rs_CO2)  = 1/MW_s(1,Carbon) *A_p*p_CO2*(r_diff*r_kin)/(r_diff+r_kin)

      else  !2020/4/17 added  
         DES_RATES(Rs_CO2)=ZERO
      endif

      write(*,*) 'DES_rates(Rs_CO2)=', DES_rates(Rs_CO2)
      write(*,*) 'Rs_CO2=', Rs_CO2
         
      
! C + H2O --> CO + H2  (kg-mole/m^3.s)
!---------------------------------------------------------------------//  
      if (DES_X_s(NP,Carbon) > c_Limiter) then    !2020/4/17 added   
         C_i=5d-12
         A_i=45.6
         E_i=4.37d7

         r_kin   = A_i *xTs * exp(-E_i/(8.3147295*1000*xTs))
         r_diff  = C_i * ((xTs+xTg)/2)**0.75/Dp0
         DES_RATES(Rs_H2O)  = 1/MW_s(1,Carbon) * A_p*p_H2O*(r_diff*r_kin)/(r_diff+r_kin)

      else  !2020/4/17 added  
         DES_RATES(Rs_H2O)  =ZERO
      endif

      write(*,*) 'DES_rates(Rs_H2O)=', DES_rates(Rs_H2O)
      write(*,*) 'Rs_H2O=', Rs_H2O
         
            
! pyrolysis, Boujjat, 2019
!---------------------------------------------------------------------//
! Biomass --> 2.77 CO + 2.71 H2 + 0.35 CO2 + 0.45 CH4 + 0.44 H2O + 0.32 Tar + 0.52 C

! meld Ku
!    chem_eq = "BiomassCore --> 0.962*CO + 0.877*H2 + 0.348*CO2 + 1.627*CH4" &
!  " + 1.421*H2O + 0.35555*Carbon"

! pyrolysis, Yin, 2020
!---------------------------------------------------------------------//
! BiomassCore --> 6.53 CO + 2.95 CO2 + 2.79 CH4 + 10.245 H2 + 6.19 C
! C18.46 H31.65 O12.43

      if (DES_X_s(NP,BiomassCore) > c_Limiter) then
         A_i=5d6     !   Ku,2015==>Prakash,2008; Also in Wang shuai,2019
         E_i=1.2d5   !   Original=1.2d8, reduced with 8314;
         DES_RATES(Pyrolysis) = A_i * exp(-E_i/(8.3147295*xTs)) * pmass(NP)*DES_X_s(NP,BiomassCore)/MW_s(1,BiomassCore) *1e-2 
         !1e-2 only for debug, please remove when it is good
      else 
         DES_RATES(Pyrolysis)=ZERO
      endif
      
      write(*,*) 'DES_rates(Pyrolysis)=', DES_rates(Pyrolysis)
      write(*,*) 'pmass(NP)=', pmass(NP)
      write(*,*) 'Pyrolysis=', Pyrolysis
      
      write(*,*) 'BiomassCore=', BiomassCore
      write(*,*) 'X_s(IJK,1,BiomassCore)=', X_s(IJK,1,BiomassCore)
      write(*,*) 'DES_X_s(NP,BiomassCore)=', DES_X_s(NP,BiomassCore)
            
   endif      !!end of if(.not. compare(EP_s(IJK,1),ZERO)) then 
      

            
      LB= NO_OF_RXNS+1     !LB=2
      UB= NO_OF_RXNS+NO_OF_DES_RXNS
      !ReactionRates(IJK,LB:UB) = ReactionRates(IJK,LB:UB) + DES_RATES
      !ReactionRates(IJK,LB:UB) =  DES_RATES
      ReactionRates(IJK,LB) =  ReactionRates(IJK,LB)     + DES_RATES(Rs_CO2)!3
      ReactionRates(IJK,LB+1) =  ReactionRates(IJK,LB+1) + DES_RATES(Rs_H2O)!4
      ReactionRates(IJK,LB+2) =  ReactionRates(IJK,LB+2) + DES_RATES(Pyrolysis)!5
            
      endif !endif of pmass>1d-10
      ENDIF !ENDIF of IF PM==1
         
      RETURN

      END SUBROUTINE USR_RATES_DES
      
