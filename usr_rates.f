      SUBROUTINE USR_RATES(IJK, RATES)

      USE compar
      USE constant
      USE energy
      USE fldvar
      USE fun_avg
      USE functions
      USE funits
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE sendrecv
      USE toleranc
      USE usr

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Gas phase concentrations (mol/m^3)
      DOUBLE PRECISION  c_H2, c_CH4, c_H2O
! Partial pressures (KPa)
      DOUBLE PRECISION  p_H2, p_H2O

! Prevent reactions when mass fraction falls below c_Limiter.
      !DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-7

!`````````````````````````````````````````````````````````````````````//
      INCLUDE 'species.inc'
!`````````````````````````````````````````````````````````````````````//


! Initialize partial pressures and concentrations.

      c_H2    = ZERO;   p_H2   = ZERO;
      c_H2O   = ZERO;   p_H2O  = ZERO;
      c_CH4   = ZERO;

! Calculate the partial pressure and molar concentration of H2
      !IF(X_g(IJK,H2) > c_Limiter) then
         c_H2 = RO_g(IJK) * X_g(IJK,H2) / Mw_g(H2)
      !ENDIF

! Calculate the molar concentration of H2O
      !IF(X_g(IJK,H2O) > c_Limiter) then
            c_H2O = RO_g(IJK)*X_g(IJK,H2O)/Mw_g(H2O)
      !ENDIF

! Calculate the molar concentration of CH4
      !IF(X_g(IJK,CH4) > c_Limiter) then
            c_CH4 = RO_g(IJK)*X_g(IJK,CH4)/Mw_g(CH4)
      !ENDIF      

!**********************************************************************!
!                         Homogeneous Reactions                        !
!**********************************************************************!
! CO + H2O <--> CO2 + H2  (kmole/m^3.s)
!---------------------------------------------------------------------//
      RATES(Rg_CO) = 0.029 * exp( 3.40d7/(8314*T_g(IJK) ) )

! CH4 + H2O --> CO + 3*H2  (kmole/m^3.s)
!---------------------------------------------------------------------//
      RATES(Rg_CH4)= 3d8 * c_CH4 * c_H2O * EXP(-1.26d8/8314/T_g(IJK))
      

! Save reaction rates for visualization

      !IF(nrr>=Rg1)  ReactionRates(IJK,Rg1)  = RATES(Rg1)
      ReactionRates(IJK,1) = RATES(Rg_CO)
      ReactionRates(IJK,2) = RATES(Rg_CH4)
      RETURN

      END SUBROUTINE USR_RATES
