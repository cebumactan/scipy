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
!           PAR(5) : EPS0
!           PAR(6) : EPS1
!           PAR(7) : X01
!           PAR(8) : 
!           PAR(9) : 
!           PAR(10): X02
!           PAR(11): 
!           PAR(12): 
!           PAR(13): X03
!           PAR(14): 
!           PAR(15): 
!           PAR(16): X11
!           PAR(17): 
!           PAR(18): 
!           PAR(19): X12
!           PAR(20): 
!           PAR(21): 
!           PAR(22): c0 coefficient 
!           PAR(23): 
!           PAR(24): 
!           PAR(25): c1 coefficient 
!           PAR(26): 

!           PAR(27): 
!           PAR(28): 
!           PAR(29): 
!           PAR(31): 
!           PAR(32): 
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
		
	D=1. + 2*A - M - N
	a=(2.+2.*A-N)/D + 2.*(1.+A)*L/D
	b=(1.+m)    /D + (1.+M+N)*L/D

        P=U(1)
	Q=U(2)
	R=U(3)

        F(1)= P * ( (R-a)/L + 2. - L*P*R - Q  )
        F(2)= Q * (  1. - L*P*R - Q  ) + b*P*R
        F(3)= R * ( (R-a)*(A-M-N)/(L*(1.+A)) + L*P*R + Q  ) / N

      END SUBROUTINE FFFF

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
! Local
      INTEGER, PARAMETER :: NDM=3
!      DOUBLE PRECISION V0(NDM),V1(NDM),G0(NDM),G1(NDM)
      DOUBLE PRECISION MAT0(NDM,NDM),MAT1(NDM,NDM)

      DOUBLE PRECISION EPS0, EPS1, A, M, N, L, D, a, b, R0, R1
      DOUBLE PRECISION mu01, mu02, mu03, mu11, mu12
      DOUBLE PRECISION X01(NDM), X01(NDM), X01(NDM), X01(NDM)
      DOUBLE PRECISION c01, c02, c03, c11, c12


        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)

	EPS0=PAR(5)
	EPS1=PAR(6)
		
	D=1. + 2*A - M - N
	a=(2.+2.*A-N)/D + 2.*(1.+A)*L/D
	b=(1.+m)    /D + (1.+M+N)*L/D

	R0 = a
	R1 = R0 - (1.+A)*L/(A-M-N)

! Linearized Matrix at M0
	MAT0(1,1) = 2.
	MAT0(1,2) = 0.
	MAT0(1,3) = 0.
	MAT0(2,1) = b*R0
	MAT0(2,2) = 1.
	MAT0(2,3) = 0.
	MAT0(3,1) = (R0*L)*R0/N
	MAT0(3,2) = R0/N
	MAT0(3,3) = (  (A-M-N)/(L*(1+A)) - N*A/( L*(1+A)*R0 )   )*R0/N

! Linearized Matrix at M1
	MAT1(1,1) = -(1+M+N)/(A-M-N)
	MAT1(1,2) = 0.
	MAT1(1,3) = 0.
	MAT1(2,1) = (b-L)*R1
	MAT1(2,2) = -1.
	MAT1(2,3) = 0.
	MAT1(3,1) = (R1*L)*R1/N
	MAT1(3,2) = R1/N
	MAT1(3,3) = (  (A-M-N)/(L*(1+A)) - N*A/( L*(1+A)*R1 )   )*R1/N

! positive eigenvalue at M0
	mu01 = 2.
	mu02 = 1.
	mu03 = -(M+N)*a/(N*L)

! negative eigenvalue at M1
	mu11 = -(1.+M+N)/(A-M-N)
	mu12 = -1.

! let AUTO compute the eigenvectors
!	FB(1)= MAT0(1,1)*X01(1) + MAT0(1,2)*X01(2) + MAT0(1,3)*X01(3)- mu01*X01(1)
!	FB(2)= MAT0(2,1)*X01(1) + MAT0(2,2)*X01(2) + MAT0(2,3)*X01(3)- mu01*X01(2)
!	FB(3)= MAT0(3,1)*X01(1) + MAT0(3,2)*X01(2) + MAT0(3,3)*X01(3)- mu01*X01(3)
!
!	FB(4)= MAT0(1,1)*X02(1) + MAT0(1,2)*X02(2) + MAT0(1,3)*X03(3)- mu01*X02(1)
!	FB(5)= MAT0(2,1)*X02(1) + MAT0(2,2)*X02(2) + MAT0(2,3)*X03(3)- mu01*X02(2)
!	FB(6)= MAT0(3,1)*X02(1) + MAT0(3,2)*X02(2) + MAT0(3,3)*X03(3)- mu01*X02(3)
!
!	FB(7)= MAT0(1,1)*X03(1) + MAT0(1,2)*X03(2) + MAT0(1,3)*X03(3)- mu01*X03(1)
!	FB(8)= MAT0(2,1)*X03(1) + MAT0(2,2)*X03(2) + MAT0(2,3)*X03(3)- mu01*X03(2)
!	FB(9)= MAT0(3,1)*X03(1) + MAT0(3,2)*X03(2) + MAT0(3,3)*X03(3)- mu01*X03(3)
!
!	FB(10)= MAT0(1,1)*X11(1) + MAT0(1,2)*X11(2) + MAT0(1,3)*X11(3)- mu01*X11(1)
!	FB(11)= MAT0(2,1)*X11(1) + MAT0(2,2)*X11(2) + MAT0(2,3)*X11(3)- mu01*X11(2)
!	FB(12)= MAT0(3,1)*X11(1) + MAT0(3,2)*X11(2) + MAT0(3,3)*X11(3)- mu01*X11(3)
!
!	FB(13)= MAT0(1,1)*X12(1) + MAT0(1,2)*X12(2) + MAT0(1,3)*X12(3)- mu01*X12(1)
!	FB(14)= MAT0(2,1)*X12(1) + MAT0(2,2)*X12(2) + MAT0(2,3)*X12(3)- mu01*X12(2)
!	FB(15)= MAT0(3,1)*X12(1) + MAT0(3,2)*X12(2) + MAT0(3,3)*X12(3)- mu01*X12(3)

! let AUTO compute the eigenvectors
	FB(1)= MAT0(1,1)*PAR(7) + MAT0(1,2)*PAR(8) + MAT0(1,3)*PAR(9)- mu01*PAR(7)
	FB(2)= MAT0(2,1)*PAR(7) + MAT0(2,2)*PAR(8) + MAT0(2,3)*PAR(9)- mu01*PAR(8)
	FB(3)= MAT0(3,1)*PAR(7) + MAT0(3,2)*PAR(8) + MAT0(3,3)*PAR(9)- mu01*PAR(9)

	FB(4)= MAT0(1,1)*PAR(10) + MAT0(1,2)*PAR(11) + MAT0(1,3)*PAR(12)- mu02*PAR(10)
	FB(5)= MAT0(2,1)*PAR(10) + MAT0(2,2)*PAR(11) + MAT0(2,3)*PAR(12)- mu02*PAR(11)
	FB(6)= MAT0(3,1)*PAR(10) + MAT0(3,2)*PAR(11) + MAT0(3,3)*PAR(12)- mu02*PAR(12)

	FB(7)= MAT0(1,1)*PAR(13) + MAT0(1,2)*PAR(14) + MAT0(1,3)*PAR(15)- mu03*PAR(13)
	FB(8)= MAT0(2,1)*PAR(13) + MAT0(2,2)*PAR(14) + MAT0(2,3)*PAR(15)- mu03*PAR(14)
	FB(9)= MAT0(3,1)*PAR(13) + MAT0(3,2)*PAR(14) + MAT0(3,3)*PAR(15)- mu03*PAR(15)

	FB(10)= MAT1(1,1)*PAR(16) + MAT1(1,2)*PAR(17) + MAT1(1,3)*PAR(18)- mu11*PAR(16)
	FB(11)= MAT1(2,1)*PAR(16) + MAT1(2,2)*PAR(17) + MAT1(2,3)*PAR(18)- mu11*PAR(17)
	FB(12)= MAT1(3,1)*PAR(16) + MAT1(3,2)*PAR(17) + MAT1(3,3)*PAR(18)- mu11*PAR(18)

	FB(13)= MAT1(1,1)*PAR(19) + MAT1(1,2)*PAR(20) + MAT1(1,3)*PAR(21)- mu12*PAR(19)
	FB(14)= MAT1(2,1)*PAR(19) + MAT1(2,2)*PAR(20) + MAT1(2,3)*PAR(21)- mu12*PAR(20)
	FB(15)= MAT1(3,1)*PAR(19) + MAT1(3,2)*PAR(20) + MAT1(3,3)*PAR(21)- mu12*PAR(21)

! unit length eigenvectors
	FB(16)= PAR(7)**2  + PAR(8)**2  + PAR(9)**2  -1
	FB(17)= PAR(10)**2 + PAR(11)**2 + PAR(12)**2 -1
	FB(18)= PAR(13)**2 + PAR(14)**2 + PAR(15)**2 -1
	FB(19)= PAR(16)**2 + PAR(17)**2 + PAR(18)**2 -1
	FB(20)= PAR(19)**2 + PAR(20)**2 + PAR(21)**2 -1

! unit length coefficients
	FB(21)= PAR(22)**2 + PAR(23)**2 + PAR(24)**2 -1
	FB(22)= PAR(25)**2 + PAR(26)**2 -1

! boundary values
	FB(23)= U0(1) - EPS0*(  PAR(22)*PAR(7) + PAR(23)*PAR(10) + PAR(24)*PAR(13)   ) 
	FB(24)= U0(2) - EPS0*(  PAR(22)*PAR(8) + PAR(23)*PAR(11) + PAR(24)*PAR(14)   ) 
	FB(25)= U0(3) - EPS0*(  PAR(22)*PAR(9) + PAR(23)*PAR(12) + PAR(24)*PAR(15)   ) 

	FB(26)= U1(1) - EPS1*(  PAR(25)*PAR(16) + PAR(26)*PAR(19)   ) 
	FB(27)= U1(2) - EPS1*(  PAR(25)*PAR(17) + PAR(26)*PAR(20)   ) 
	FB(28)= U1(3) - EPS1*(  PAR(25)*PAR(18) + PAR(26)*PAR(21)   ) 


      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS





