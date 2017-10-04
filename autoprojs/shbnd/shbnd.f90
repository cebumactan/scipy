!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   shbnd :    Heteroclinic orbits : A saddle-saddle connection when alpha=0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameter assignment:
!
!           PAR(1) : A = alpha
!           PAR(2) : M = m
!           PAR(3) : N = n
!           PAR(4) : L = lambda
!           PAR(5) :
!           PAR(6) : 
!           PAR(7) : 
!           PAR(8) : 
!           PAR(9) :
!           PAR(10): 
!           PAR(11): 
!           PAR(12): 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION PERIOD, DUMMY(1)

	CALL FFFF(3,U,ICP,PAR,0,F,DUMMY)

	PERIOD=PAR(11)
	F(1)=PERIOD*F(1)
	F(2)=PERIOD*F(2)
	F(3)=PERIOD*F(3)

      END SUBROUTINE FUNC

      SUBROUTINE FFFF(NDM,U,ICP,PAR,IJAC,F,DFDU)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDM,NDM)

      DOUBLE PRECISION A, M, N, L, D, a, b, P, Q, R

        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)
		
	D=1 + 2*A - M - N
	a=(2+2*A-N)/D + 2*(1+A)*L/D
	b=(1+m)    /D + (1+M+N)*L/D

        P=U(1)
	Q=U(2)
	R=U(3)

        F(1)= P * ( (R-a)/L + 2 - L*P*R - Q  )
        F(2)= Q * (  1 - L*P*R - Q  ) + b*P*R
        F(3)= R * ( (R-a)*(A-M-N)/(L*(1+A)) + L*P*R + Q  ) / N

      END SUBROUTINE FFFF
