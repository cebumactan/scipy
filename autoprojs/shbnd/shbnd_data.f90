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
!           PAR(27): PERIOD
!           PAR(28): mu01
!           PAR(29): mu02
!           PAR(30): mu03
!           PAR(31): mu11
!           PAR(32): mu12
!	    PAR(33): K    My Convenient Constants From Here
!	    PAR(34): E0
!	    PAR(35): F0
!	    PAR(36): C1
!	    PAR(37): Q00
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

	PERIOD=PAR(27)
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

      DOUBLE PRECISION A, M, N, L, D, aa, bb, P, Q, R

        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+M)    /D + (1.D0+M+N)*L/D

        P=U(1)
	Q=U(2)
	R=U(3)

        F(1)= P * ( (R-aa)/L + 2.D0 - L*P*R - Q  )
        F(2)= Q * (  1.D0 - L*P*R - Q  ) + bb*P*R
        F(3)= R * ( (R-aa)*(A-M-N)/(L*(1.D0+A)) + L*P*R + Q  ) / N

      END SUBROUTINE FFFF



      SUBROUTINE FAC(K,VAL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K
      DOUBLE PRECISION, INTENT(OUT) :: VAL
      INTEGER j
        VAL=1.D0
	DO j=1,K
	    VAL=VAL*j
	END DO
      END SUBROUTINE FAC

      SUBROUTINE COMB(K,J,VAL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K,J
      DOUBLE PRECISION, INTENT(OUT) :: VAL
      DOUBLE PRECISION DUMMY
	CALL FAC(K,DUMMY)
	VAL = DUMMY
	CALL FAC(J,DUMMY)
	VAL = VAL/ DUMMY
	CALL FAC(K-J,DUMMY)
	VAL = VAL/ DUMMY
      END SUBROUTINE COMB

! EXPLICIT SOLUTION FOR r(x)
      SUBROUTINE RRRR(X,R,PAR)
!     ---------- ----

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: X, PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: R
      DOUBLE PRECISION A, M, N, L, D, aa, bb, PERIOD
      DOUBLE PRECISION R0, R1, Q00, R00, EPS0, NUMERATOR, DENOM, Z, DUMMY,E0, F0, C1, XMAX
      INTEGER K,j

        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)
	K=PAR(33)
	C1=PAR(36)

	PERIOD = PAR(27)
	XMAX = 0.5*PERIOD
		
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+m)    /D + (1.D0+M+N)*L/D

! Equilibrium points
	R0 = aa
	R1 = R0 - (1.D0+A)*L/(A-M-N)

! Some Constants
	Z = -aa*(M+N)/L
	E0 = C1/(C1 + DEXP(XMAX))
	F0 = 1.0D0/(C1 + DEXP(XMAX))

	NUMERATOR = aa*( (E0+F0*DEXP(X))**K )

	DENOM = 0.D0
	DO j = 0, K
		CALL COMB(K,j,DUMMY)
		DENOM = DENOM +  ( DUMMY* ( E0**(K-j) ) * ( F0**j ) * DEXP(j*X) ) / (-K*Z+j)
	END DO
	DENOM = DENOM*(-K*Z)


	R = NUMERATOR / DENOM

      END SUBROUTINE RRRR



! EXPLICIT SOLUTION
! T runs in [0,1]			: SCALED VAR
! X runs in [-0.5*PERIOD, 0.5*PERIOD]	: ORIGINAL VAR
! i.e., X = PERIOD*(T-0.5)

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      INTEGER, PARAMETER :: NDM=3
      DOUBLE PRECISION PERIOD, X, EPS0, EPS1
      DOUBLE PRECISION A, M, N, L, LMAX, D, aa, bb, P, Q, R
      DOUBLE PRECISION R0, R1, S0, S1, mu01, mu02, mu03, mu11, mu12, DUMMY
      DOUBLE PRECISION X01(NDM), X02(NDM), X03(NDM), X11(NDM), X12(NDM), VEC(NDM) 
      DOUBLE PRECISION C1, P00, Q00, R00, P11, Q11, R11, XMAX, XMIN, Z, E0, F0, CC0, CC1, CC2
      INTEGER K,j


! INITIALIZE 

! Try Interval Length 20
        PERIOD=20.D0

        A=0.D0
        M=-0.5D0
	K=10
	N=1.D0/K
	LMAX = 2*(A-M-N)*(1+M)/((1+M+N)**2)
	L=0.5*LMAX
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+M)    /D + (1.D0+M+N)*L/D

! Equilibrium points
	R0 = aa
	R1 = R0 - (1.D0+A)*L/(A-M-N)
	S0 = (1+M+N)/(1+A) - N/(R0*(1+A))
	S1 = (1+M+N)/(1+A) - N/(R1*(1+A))

! positive eigenvaalue at M0
	mu01 = 2.D0
	mu02 = 1.D0
	mu03 = -(M+N)*aa/(N*L)

! negative eigenvalue at M1
	mu11 = -(1.D0+M+N)/(A-M-N)
	mu12 = -1.D0


! Provide the eigenvectors with unit length
	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu01/S0 ) - (N/R0)*( 1.D0/L + mu01 )*( R0/L + mu01/S0 )
	 X01(1)=1.D0
	 X01(2)=bb*R0
	 X01(3)=-(L+bb)*R0* ( (1.D0+A)*R0/L + mu01/S0  ) / DUMMY
         DUMMY = SQRT( X01(1)**2 + X01(2)**2 + X01(3)**2 )
	 X01(1)=X01(1)/DUMMY
	 X01(2)=X01(2)/DUMMY
	 X01(3)=X01(3)/DUMMY

	 WRITE(*,*) X01(1), X01(2), X01(3)

	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu02/S0 ) - (N/R0)*( 1.D0/L + mu02 )*( R0/L + mu02/S0 )
	 X02(1)=0.D0
	 X02(2)=1.D0
	 X02(3)=-( (1.D0+A)*R0/L + mu02/S0 )/DUMMY
         DUMMY = SQRT( X02(1)**2 + X02(2)**2 + X02(3)**2 )
	 X02(1)=X02(1)/DUMMY
	 X02(2)=X02(2)/DUMMY
	 X02(3)=X02(3)/DUMMY

	 X03(1)=0.D0
	 X03(2)=0.D0
	 X03(3)=1.D0

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu11/s1 ) - (N/R1)*( 1.D0/L + mu11 )*( R1/L + mu11/S1 )
	 X11(1)=1.D0
	 X11(2)=(bb-L)*R1/(1.D0+mu11)
	 X11(3)=-( L*R1 + X11(2) ) * ( (1.D0+A)*R1/L + mu11/S1 ) / DUMMY
         DUMMY = SQRT( X11(1)**2 + X11(2)**2 + X11(3)**2 )
	 X11(1)=X11(1)/DUMMY
	 X11(2)=X11(2)/DUMMY
	 X11(3)=X11(3)/DUMMY

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu12/s1 ) - (N/R1)*( 1.D0/L + mu12 )*( R1/L + mu12/S1 )
	 X12(1)=0.D0
	 X12(2)=1.D0
	 X12(3)=- ( (1.D0+A)*R1/L + mu12/S1 ) / DUMMY
         DUMMY = SQRT( X12(1)**2 + X12(2)**2 + X12(3)**2 )
	 X12(1)=X12(1)/DUMMY
	 X12(2)=X12(2)/DUMMY
	 X12(3)=X12(3)/DUMMY


	XMAX = 0.5D0*PERIOD
	XMIN = -XMAX

	C1 = 1.D0
	P00 = 0.D0
	Q00 = 1.D0/( 1.D0 + C1*DEXP(-(XMIN)) )
	P11 = 0.D0
	Q11 = 1.D0/( 1.D0 + C1*DEXP(-(XMAX)) )
	E0 = C1
	F0 = 1.D0

!	WRITE(*,*) mu01, mu02, mu03, mu11, mu12

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

!           PAR(27): PERIOD
!           PAR(28): mu01
!           PAR(29): mu02
!           PAR(30): mu03
!           PAR(31): mu11
!           PAR(32): mu12

!	    PAR(33): K    My Convenient Constants From Here
!	    PAR(34): E0
!	    PAR(35): F0
!	    PAR(36): C1
!	    PAR(37): Q00

         PAR(1) = A
         PAR(2) = M
         PAR(3) = N
         PAR(4) = L
         PAR(5) = 1.0488740383E-08
 	 PAR(6) = 7.0036896621D-05

         PAR(7) = X01(1) 
         PAR(8) = X01(2)
         PAR(9) = X01(3)

         PAR(10)= X02(1)
         PAR(11)= X02(2)
         PAR(12)= X02(3)

         PAR(13)= X03(1)
         PAR(14)= X03(2)
         PAR(15)= X03(3)

         PAR(16)= X11(1)
         PAR(17)= X11(2)
         PAR(18)= X11(3)

         PAR(19)= X12(1)
         PAR(20)= X12(2)
         PAR(21)= X12(3)

         PAR(22)= 9.8572084967D-01
         PAR(23)= 1.6838766735D-01
         PAR(24)= 1.5512729904D-06 

         PAR(25)= 5.4228550939D-03
         PAR(26)= -9.9998529621D-01

         PAR(27)= PERIOD
         PAR(28)= mu01
         PAR(29)= mu02
         PAR(30)= mu03
         PAR(31)= mu11
         PAR(32)= mu12
	 PAR(33)= K    

	 PAR(36)= C1



      END SUBROUTINE STPNT



! Constraints
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
! Local
      INTEGER, PARAMETER :: NDM=3
      DOUBLE PRECISION MAT0(NDM,NDM),MAT1(NDM,NDM)

      DOUBLE PRECISION EPS0, EPS1, A, M, N, L, D, aa, bb, R0, R1
      DOUBLE PRECISION mu01, mu02, mu03, mu11, mu12
      DOUBLE PRECISION X01(NDM), X02(NDM), X03(NDM), X11(NDM), X12(NDM)

        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)

	EPS0=PAR(5)
	EPS1=PAR(6)
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+M)    /D + (1.D0+M+N)*L/D

	R0 = aa
	R1 = R0 - (1.D0+A)*L/(A-M-N)

! unit length coefficients
	FB(1)= PAR(22)**2 + PAR(23)**2 + PAR(24)**2 -1
	FB(2)= PAR(25)**2 + PAR(26)**2 -1

! boundary values
	FB(3)= U0(1) - ( 0.0D0 + EPS0*(  PAR(22)*PAR(7) + PAR(23)*PAR(10) + PAR(24)*PAR(13)   ) )
	FB(4)= U0(2) - ( 0.0D0 + EPS0*(  PAR(22)*PAR(8) + PAR(23)*PAR(11) + PAR(24)*PAR(14)   ) )
	FB(5)= U0(3) - ( R0    + EPS0*(  PAR(22)*PAR(9) + PAR(23)*PAR(12) + PAR(24)*PAR(15)   ) )

	FB(6)= U1(1) - ( 0.0D0 + EPS1*(  PAR(25)*PAR(16) + PAR(26)*PAR(19)   ) )
	FB(7)= U1(2) - ( 1.0D0 + EPS1*(  PAR(25)*PAR(17) + PAR(26)*PAR(20)   ) )
	FB(8)= U1(3) - ( R1    + EPS1*(  PAR(25)*PAR(18) + PAR(26)*PAR(21)   ) )


!
! Provide the eigenvectors with unit length
!	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu01/S0 ) - (N/R0)*( 1.D0/L + mu01 )*( R0/L + mu01/S0 )
!	 X01(1)=1.D0
!	 X01(2)=bb*R0
!	 X01(3)=-(L+bb)*R0* ( (1.D0+A)*R0/L + mu01/S0  ) / DUMMY
!         DUMMY = SQRT( X01(1)**2 + X01(2)**2 + X01(3)**2 )
!	 X01(1)=X01(1)/DUMMY
!	 X01(2)=X01(2)/DUMMY
!	 X01(3)=X01(3)/DUMMY
!
!	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu02/S0 ) - (N/R0)*( 1.D0/L + mu02 )*( R0/L + mu02/S0 )
!	 X02(1)=0.D0
!	 X02(2)=1.D0
!	 X02(3)=-( (1.D0+A)*R0/L + mu02/S0 )/DUMMY
!         DUMMY = SQRT( X02(1)**2 + X02(2)**2 + X02(3)**2 )
!	 X02(1)=X02(1)/DUMMY
!	 X02(2)=X02(2)/DUMMY
!	 X02(3)=X02(3)/DUMMY
!
!	 X03(1)=0.D0
!	 X03(2)=0.D0
!	 X03(3)=1.D0
!
!	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu11/s1 ) - (N/R1)*( 1.D0/L + mu11 )*( R1/L + mu11/S1 )
!	 X11(1)=1.D0
!	 X11(2)=(B-L)*R1/(1.D0+mu11)
!	 X11(3)=-( L*R1 + X11(2) ) * ( (1.D0+A)*R1/L + mu11/S1 ) / DUMMY
!        DUMMY = SQRT( X11(1)**2 + X11(2)**2 + X11(3)**2 )
!	 X11(1)=X11(1)/DUMMY
!	 X11(2)=X11(2)/DUMMY
!	 X11(3)=X11(3)/DUMMY
!
!	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu12/s1 ) - (N/R1)*( 1.D0/L + mu12 )*( R1/L + mu12/S1 )
!	 X12(1)=0.D0
!	 X12(2)=1.D0
!	 X12(3)=- ( (1.D0+A)*R1/L + mu12/S1 ) / DUMMY
!        DUMMY = SQRT( X12(1)**2 + X12(2)**2 + X12(3)**2 )
!	 X12(1)=X12(1)/DUMMY
!	 X12(2)=X12(2)/DUMMY
!	 X12(3)=X12(3)/DUMMY
!
!
! Linearized Matrix at M0
!	MAT0(1,1) = 2.D0
!	MAT0(1,2) = 0.D0
!	MAT0(1,3) = 0.D0
!	MAT0(2,1) = bb*R0
!	MAT0(2,2) = 1.D0
!	MAT0(2,3) = 0.D0
!	MAT0(3,1) = (R0*L)*R0/N
!	MAT0(3,2) = R0/N
!	MAT0(3,3) = (  (A-M-N)/(L*(1.D0+A)) - N*A/( L*(1.D0+A)*R0 )   )*R0/N
!
! Linearized Matrix at M1
!	MAT1(1,1) = -(1.D0+M+N)/(A-M-N)
!	MAT1(1,2) = 0.D0
!	MAT1(1,3) = 0.D0
!	MAT1(2,1) = (bb-L)*R1
!	MAT1(2,2) = -1.D0
!	MAT1(2,3) = 0.D0
!	MAT1(3,1) = (R1*L)*R1/N
!	MAT1(3,2) = R1/N
!	MAT1(3,3) = (  (A-M-N)/(L*(1.D0+A)) - N*A/( L*(1.D0+A)*R1 )   )*R1/N
!
! positive eigenvalue at M0
!	mu01 = 2.D0
!	mu02 = 1.D0
!	mu03 = -(M+N)*aa/(N*L)
!
! negative eigenvalue at M1
!	mu11 = -(1.D0+M+N)/(A-M-N)
!	mu12 = -1.D0
!
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
!	FB(1)= MAT0(1,1)*PAR(7) + MAT0(1,2)*PAR(8) + MAT0(1,3)*PAR(9)- mu01*PAR(7)
!	FB(2)= MAT0(2,1)*PAR(7) + MAT0(2,2)*PAR(8) + MAT0(2,3)*PAR(9)- mu01*PAR(8)
!	FB(3)= MAT0(3,1)*PAR(7) + MAT0(3,2)*PAR(8) + MAT0(3,3)*PAR(9)- mu01*PAR(9)
!
!	FB(4)= MAT0(1,1)*PAR(10) + MAT0(1,2)*PAR(11) + MAT0(1,3)*PAR(12)- mu02*PAR(10)
!	FB(5)= MAT0(2,1)*PAR(10) + MAT0(2,2)*PAR(11) + MAT0(2,3)*PAR(12)- mu02*PAR(11)
!	FB(6)= MAT0(3,1)*PAR(10) + MAT0(3,2)*PAR(11) + MAT0(3,3)*PAR(12)- mu02*PAR(12)
!
!	FB(7)= MAT0(1,1)*PAR(13) + MAT0(1,2)*PAR(14) + MAT0(1,3)*PAR(15)- mu03*PAR(13)
!	FB(8)= MAT0(2,1)*PAR(13) + MAT0(2,2)*PAR(14) + MAT0(2,3)*PAR(15)- mu03*PAR(14)
!	FB(9)= MAT0(3,1)*PAR(13) + MAT0(3,2)*PAR(14) + MAT0(3,3)*PAR(15)- mu03*PAR(15)
!
!	FB(10)= MAT1(1,1)*PAR(16) + MAT1(1,2)*PAR(17) + MAT1(1,3)*PAR(18)- mu11*PAR(16)
!	FB(11)= MAT1(2,1)*PAR(16) + MAT1(2,2)*PAR(17) + MAT1(2,3)*PAR(18)- mu11*PAR(17)
!	FB(12)= MAT1(3,1)*PAR(16) + MAT1(3,2)*PAR(17) + MAT1(3,3)*PAR(18)- mu11*PAR(18)
!
!	FB(13)= MAT1(1,1)*PAR(19) + MAT1(1,2)*PAR(20) + MAT1(1,3)*PAR(21)- mu12*PAR(19)
!	FB(14)= MAT1(2,1)*PAR(19) + MAT1(2,2)*PAR(20) + MAT1(2,3)*PAR(21)- mu12*PAR(20)
!	FB(15)= MAT1(3,1)*PAR(19) + MAT1(3,2)*PAR(20) + MAT1(3,3)*PAR(21)- mu12*PAR(21)
! unit length eigenvectors
!	FB(16)= PAR(7)**2  + PAR(8)**2  + PAR(9)**2  -1
!	FB(17)= PAR(10)**2 + PAR(11)**2 + PAR(12)**2 -1
!	FB(18)= PAR(13)**2 + PAR(14)**2 + PAR(15)**2 -1
!	FB(19)= PAR(16)**2 + PAR(17)**2 + PAR(18)**2 -1
!	FB(20)= PAR(19)**2 + PAR(20)**2 + PAR(21)**2 -1
! use explicit formulas for eigenvectors
!	FB(9)=PAR(7)-X01(1)
!	FB(10)=PAR(8)-X01(2)
!	FB(11)=PAR(9)-X01(3)
!
!	FB(12)=PAR(10)-X02(1)
!	FB(13)=PAR(11)-X02(2)
!	FB(14)=PAR(12)-X02(3)
!
!	FB(15)=PAR(13)-X03(1)
!	FB(16)=PAR(14)-X03(2)
!	FB(17)=PAR(15)-X03(3)
!
!	FB(18)=PAR(16)-X11(1)
!	FB(19)=PAR(17)-X11(2)
!	FB(20)=PAR(18)-X11(3)
!
!	FB(21)=PAR(19)-X12(1)
!	FB(22)=PAR(20)-X12(2)
!	FB(23)=PAR(21)-X12(3)
!
! use explicit formulas for eigenvalues
!	FB(24)=PAR(28)-mu01
!	FB(25)=PAR(29)-mu02
!	FB(26)=PAR(30)-mu03
!	FB(27)=PAR(31)-mu11
!	FB(28)=PAR(32)-mu12


      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS





