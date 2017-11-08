!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   shbnd :    Heteroclinic orbits : A saddle-saddle connection 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameter assignment:
!
!           PAR(1) : A = alpha
!           PAR(2) : M = m
!           PAR(3) : N = n
!           PAR(4) : L = lambda
!           PAR(5) : K integer = 1/N
!           PAR(6) : C1 
!           PAR(7) : EPS0
!           PAR(8) : EPS1
!           PAR(9) : c01
!           PAR(10): c02
!           PAR(11): c03
!           PAR(12): c11
!           PAR(13): c12
!           PAR(14): c14
!           PAR(15): PERIOD
!           PAR(16): mu01
!           PAR(17): mu02
!           PAR(18): mu03
!           PAR(19): mu11
!           PAR(20): mu12
!           PAR(21): mu14

!           PAR(22): X01
!           PAR(23): 
!           PAR(24): 
!           PAR(25): 

!           PAR(26): X02
!           PAR(27): 
!           PAR(28): 
!           PAR(29): 

!           PAR(30): X03 
!           PAR(31): 
!           PAR(32): 
!	    PAR(33): 

!	    PAR(34): X11
!	    PAR(35): 
!	    PAR(36): 
!	    PAR(37):
 
!           PAR(38): X12
!	    PAR(39): 
!	    PAR(40): 
!	    PAR(41):
 
!	    PAR(42): X14
!	    PAR(43): 
!           PAR(44): 
!	    PAR(45): 

!	    PAR(46): Extra
!	    PAR(47): 
!	    PAR(48): 
!	    PAR(49): 
!	    PAR(50):
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

	CALL FFFF(4,U,ICP,PAR,0,F,DUMMY)

	PERIOD=PAR(15)
	F(1)=PERIOD*F(1)
	F(2)=PERIOD*F(2)
	F(3)=PERIOD*F(3)
	F(4)=PERIOD*F(4)

      END SUBROUTINE FUNC

      SUBROUTINE FFFF(NDM,U,ICP,PAR,IJAC,F,DFDU)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDM,NDM)

      DOUBLE PRECISION A, M, N, L, D, aa, bb, P, Q, R, S

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
	S=U(4)

        F(1)= P * ( (R-aa)/L + 2.D0 - L*P*R - Q  )
        F(2)= Q * (  1.D0 - L*P*R - Q  ) + bb*P*R
        F(3)= R * ( (R-aa)*(A-M-N)/(L*(1.D0+A)) + L*P*R + Q +   ( R*( S-(1.D0+M+N)/(1.D0+A) ) + N/(1.D0+A) )*A/L  ) / N
	F(4)= S * ( (R-aa)*(A-M-N)/(L*(1.D0+A)) + L*P*R + Q -   ( R*( S-(1.D0+M+N)/(1.D0+A) ) + N/(1.D0+A) )  /L  )

      END SUBROUTINE FFFF


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

      INTEGER, PARAMETER :: NDM=4
      DOUBLE PRECISION PERIOD, X, EPS0, EPS1
      DOUBLE PRECISION A, M, N, L, LMAX, D, aa, bb, P, Q, R, S
      DOUBLE PRECISION R0, R1, S0, S1, mu01, mu02, mu03, mu04, mu11, mu12, mu13, mu14, DUMMY
      DOUBLE PRECISION X01(NDM), X02(NDM), X03(NDM), X04(NDM), X11(NDM), X12(NDM), X13(NDM), X14(NDM), VEC(NDM) 
      DOUBLE PRECISION C1, P00, Q00, R00, P11, Q11, R11, XMAX, XMIN, Z, E0, F0, CC0, CC1, CC2, QA, QB, QC
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

	C1 = 1.D0
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+M)    /D + (1.D0+M+N)*L/D

! Equilibrium points
	R0 = aa
	R1 = R0 - (1.D0+A)*L/(A-M-N)
	S0 = (1.D0+M+N)/(1.D0+A) - N/(R0*(1.D0+A))
	S1 = (1.D0+M+N)/(1.D0+A) - N/(R1*(1.D0+A))

! positive eigenvaalue at M0
	QA = 1.D0
	QB = -( (1.D0-S0)/L - N/(L*R0) )*R0/N + S0*R0/L
	QC = -( S0*R0*R0/(L*N) )*( (1.D0-S0)/L - N/(L*R0) ) - (S0*R0/N)*((1.D0-S0)/L)*(A*R0/L)

	mu01 = 2.D0
	mu02 = 1.D0
	mu03 = (  -QB + DSQRT(QB**2-4*QA*QC)  )/(2*QA)
	mu04 = (  -QB - DSQRT(QB**2-4*QA*QC)  )/(2*QA)

! negative eigenvalue at M1
	QA = 1.D0
	QB = -( (1.D0-S1)/L - N/(L*R1) )*R1/N + S1*R1/L
	QC = -( S1*R1*R1/(L*N) )*( (1.D0-S1)/L - N/(L*R1) ) - (S1*R1/N)*((1.D0-S1)/L)*(A*R1/L)

	mu11 = -(1.D0+M+N)/(A-M-N)
	mu12 = -1.D0
	mu13 = (  -QB + DSQRT(QB**2-4*QA*QC)  )/(2*QA)
	mu14 = (  -QB - DSQRT(QB**2-4*QA*QC)  )/(2*QA)

! Provide the eigenvectors with unit length
	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu01/S0 ) - (N/R0)*( 1.D0/L + mu01 )*( R0/L + mu01/S0 )
	 X01(1)=1.D0
	 X01(2)=bb*R0
	 X01(3)=-(L+bb)*R0* ( (1.D0+A)*R0/L + mu01/S0  ) / DUMMY
	 X01(4)=-(L+bb)*R0* ( N*(1.D0/L + mu01) / R0   ) / DUMMY
         DUMMY = DSQRT( X01(1)**2 + X01(2)**2 + X01(3)**2 + X01(4)**2)
	 X01(1)=X01(1)/DUMMY
	 X01(2)=X01(2)/DUMMY
	 X01(3)=X01(3)/DUMMY
	 X01(4)=X01(4)/DUMMY


	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu02/S0 ) - (N/R0)*( 1.D0/L + mu02 )*( R0/L + mu02/S0 )
	 X02(1)=0.D0
	 X02(2)=1.D0
	 X02(3)=-( (1.D0+A)*R0/L + mu02/S0 )/DUMMY
	 X02(4)=-( N*( 1.D0/L + mu02 )/ R0 )/DUMMY
         DUMMY = DSQRT( X02(1)**2 + X02(2)**2 + X02(3)**2 + X02(4)**2)
	 X02(1)=X02(1)/DUMMY
	 X02(2)=X02(2)/DUMMY
	 X02(3)=X02(3)/DUMMY
	 X02(4)=X02(4)/DUMMY

	 X03(1)=0.D0
	 X03(2)=0.D0
	 X03(3)=1.D0
	 X03(4)=N*( (1.D0-S0)/L  ) / ( N*R0/L + N*mu03/S0 )
         DUMMY = DSQRT( X03(1)**2 + X03(2)**2 + X03(3)**2 + X03(4)**2)
	 X03(1)=X03(1)/DUMMY
	 X03(2)=X03(2)/DUMMY
	 X03(3)=X03(3)/DUMMY
	 X03(4)=X03(4)/DUMMY

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu11/s1 ) - (N/R1)*( 1.D0/L + mu11 )*( R1/L + mu11/S1 )
	 X11(1)=1.D0
	 X11(2)=(bb-L)*R1/(1.D0+mu11)
	 X11(3)=-( L*R1 + X11(2) ) * ( (1.D0+A)*R1/L + mu11/S1 ) / DUMMY
	 X11(4)=-( L*R1 + X11(2) ) * ( N*( 1.D0/L + mu11)/R1   ) / DUMMY
         DUMMY = DSQRT( X11(1)**2 + X11(2)**2 + X11(3)**2  + X11(4)**2)
	 X11(1)=X11(1)/DUMMY
	 X11(2)=X11(2)/DUMMY
	 X11(3)=X11(3)/DUMMY
	 X11(4)=X11(4)/DUMMY

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu12/s1 ) - (N/R1)*( 1.D0/L + mu12 )*( R1/L + mu12/S1 )
	 X12(1)=0.D0
	 X12(2)=1.D0
	 X12(3)=- ( (1.D0+A)*R1/L + mu12/S1 ) / DUMMY
	 X12(4)=- ( N*( 1.D0/L + mu12 )/R1  ) /DUMMY
         DUMMY = DSQRT( X12(1)**2 + X12(2)**2 + X12(3)**2 + X12(4)**2)
	 X12(1)=X12(1)/DUMMY
	 X12(2)=X12(2)/DUMMY
	 X12(3)=X12(3)/DUMMY
	 X12(4)=X12(4)/DUMMY


	 X14(1)=0.D0
	 X14(2)=0.D0
	 X14(3)= ( R1/L + mu14/S1 ) / ( (1.D0-S1)/L )
	 X14(4)=1.D0
         DUMMY = DSQRT( X14(1)**2 + X14(2)**2 + X14(3)**2 + X14(4)**2)
	 X14(1)=X14(1)/DUMMY
	 X14(2)=X14(2)/DUMMY
	 X14(3)=X14(3)/DUMMY
	 X14(4)=X14(4)/DUMMY

	 X13(1)=0.D0
	 X13(2)=0.D0
	 X13(3)=1.D0
	 X13(4)=N* ((1.D0-S1)/L) / ( N*R1/L + N*mu13/S1)
         DUMMY = DSQRT( X13(1)**2 + X13(2)**2 + X13(3)**2 + X13(4)**2)
	 X13(1)=X13(1)/DUMMY
	 X13(2)=X13(2)/DUMMY
	 X13(3)=X13(3)/DUMMY
	 X13(4)=X13(4)/DUMMY

	 X04(1)=0.D0
	 X04(2)=0.D0
	 X04(3)= ( R0/L + mu04/S0 ) / ( (1.D0-S0)/L )
	 X04(4)=1.D0
         DUMMY = DSQRT( X04(1)**2 + X04(2)**2 + X04(3)**2 + X04(4)**2)
	 X04(1)=X04(1)/DUMMY
	 X04(2)=X04(2)/DUMMY
	 X04(3)=X04(3)/DUMMY
	 X04(4)=X04(4)/DUMMY



         PAR(1) = A
         PAR(2) = M
         PAR(3) = N
         PAR(4) = L
	 PAR(5) = K    
	 PAR(6)= C1

         PAR(7) = 8.1282436837416404E-05
 	 PAR(8) = 0.00027899278715616039

         PAR(9)= 0.E0
!         PAR(9)= 1.0D-1
         PAR(10)= 1.E0
         PAR(11)= 8.64794907E-07 

         PAR(12)= 0.E0
         PAR(13)= -2.52412665E-01
	 PAR(14)= 9.67619681E-01


         PAR(15)= PERIOD
         PAR(16)= mu01
         PAR(17)= mu02
         PAR(18)= mu03
         PAR(19)= mu11
         PAR(20)= mu12
	 PAR(21)= mu14

         PAR(22)= X01(1) 
         PAR(23)= X01(2)
         PAR(24)= X01(3)
         PAR(25)= X01(4)

         PAR(26)= X02(1)
         PAR(27)= X02(2)
         PAR(28)= X02(3)
         PAR(29)= X02(4)

         PAR(30)= X03(1)
         PAR(31)= X03(2)
         PAR(32)= X03(3)
         PAR(33)= X03(4)

         PAR(34)= X11(1)
         PAR(35)= X11(2)
         PAR(36)= X11(3)
         PAR(37)= X11(4)

         PAR(38)= X12(1)
         PAR(39)= X12(2)
         PAR(40)= X12(3)
         PAR(41)= X12(4)

         PAR(42)= X14(1)
         PAR(43)= X14(2)
         PAR(44)= X14(3)
         PAR(45)= X14(4)

	 PAR(46)= X13(1)
	 PAR(47)= X13(2)
	 PAR(48)= X13(3)
	 PAR(49)= X13(4)
	

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
      INTEGER, PARAMETER :: NDM=4
      DOUBLE PRECISION MAT0(NDM,NDM),MAT1(NDM,NDM)

      DOUBLE PRECISION EPS0, EPS1, A, M, N, L, D, aa, bb, R0, R1, S0, S1
      DOUBLE PRECISION mu01, mu02, mu03, mu04, mu11, mu12, mu13, mu14
      DOUBLE PRECISION X01(NDM), X02(NDM), X03(NDM), X04(NDM), X11(NDM), X12(NDM), X13(NDM), X14(NDM), QA, QB, QC, DUMMY

        A=PAR(1)
        M=PAR(2)
	N=PAR(3)
	L=PAR(4)

	EPS0=PAR(7)
	EPS1=PAR(8)
		
	D=1.D0 + 2.D0*A - M - N
	aa=(2.D0+2.D0*A-N)/D + 2.D0*(1.D0+A)*L/D
	bb=(1.D0+M)    /D + (1.D0+M+N)*L/D


! Equilibrium points
	R0 = aa
	R1 = R0 - (1.D0+A)*L/(A-M-N)
	S0 = (1.D0+M+N)/(1.D0+A) - N/(R0*(1.D0+A))
	S1 = (1.D0+M+N)/(1.D0+A) - N/(R1*(1.D0+A))

! positive eigenvaalue at M0
	QA = 1.D0
	QB = -( (1.D0-S0)/L - N/(L*R0) )*R0/N + S0*R0/L
	QC = -( S0*R0*R0/(L*N) )*( (1.D0-S0)/L - N/(L*R0) ) - (S0*R0/N)*((1.D0-S0)/L)*(A*R0/L)

	mu01 = 2.D0
	mu02 = 1.D0
	mu03 = (  -QB + DSQRT(QB**2-4*QA*QC)  )/(2*QA)
	mu04 = (  -QB - DSQRT(QB**2-4*QA*QC)  )/(2*QA)

! negative eigenvalue at M1
	QA = 1.D0
	QB = -( (1.D0-S1)/L - N/(L*R1) )*R1/N + S1*R1/L
	QC = -( S1*R1*R1/(L*N) )*( (1.D0-S1)/L - N/(L*R1) ) - (S1*R1/N)*((1.D0-S1)/L)*(A*R1/L)

	mu11 = -(1.D0+M+N)/(A-M-N)
	mu12 = -1.D0
	mu13 = (  -QB + DSQRT(QB**2-4*QA*QC)  )/(2*QA)
	mu14 = (  -QB - DSQRT(QB**2-4*QA*QC)  )/(2*QA)

! Provide the eigenvectors with unit length
	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu01/S0 ) - (N/R0)*( 1.D0/L + mu01 )*( R0/L + mu01/S0 )
	 X01(1)=1.D0
	 X01(2)=bb*R0
	 X01(3)=-(L+bb)*R0* ( (1.D0+A)*R0/L + mu01/S0  ) / DUMMY
	 X01(4)=-(L+bb)*R0* ( N*(1.D0/L + mu01) / R0   ) / DUMMY
         DUMMY = DSQRT( X01(1)**2 + X01(2)**2 + X01(3)**2 + X01(4)**2)
	 X01(1)=X01(1)/DUMMY
	 X01(2)=X01(2)/DUMMY
	 X01(3)=X01(3)/DUMMY
	 X01(4)=X01(4)/DUMMY


	 DUMMY = ( (1.D0-S0)/L ) * ( (1.D0+A)*R0/L + mu02/S0 ) - (N/R0)*( 1.D0/L + mu02 )*( R0/L + mu02/S0 )
	 X02(1)=0.D0
	 X02(2)=1.D0
	 X02(3)=-( (1.D0+A)*R0/L + mu02/S0 )/DUMMY
	 X02(4)=-( N*( 1.D0/L + mu02 )/ R0 )/DUMMY
         DUMMY = DSQRT( X02(1)**2 + X02(2)**2 + X02(3)**2 + X02(4)**2)
	 X02(1)=X02(1)/DUMMY
	 X02(2)=X02(2)/DUMMY
	 X02(3)=X02(3)/DUMMY
	 X02(4)=X02(4)/DUMMY

	 X03(1)=0.D0
	 X03(2)=0.D0
	 X03(3)=1.D0
	 X03(4)=N*( (1.D0-S0)/L  ) / ( N*R0/L + N*mu03/S0 )
         DUMMY = DSQRT( X03(1)**2 + X03(2)**2 + X03(3)**2 + X03(4)**2)
	 X03(1)=X03(1)/DUMMY
	 X03(2)=X03(2)/DUMMY
	 X03(3)=X03(3)/DUMMY
	 X03(4)=X03(4)/DUMMY

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu11/s1 ) - (N/R1)*( 1.D0/L + mu11 )*( R1/L + mu11/S1 )
	 X11(1)=1.D0
	 X11(2)=(bb-L)*R1/(1.D0+mu11)
	 X11(3)=-( L*R1 + X11(2) ) * ( (1.D0+A)*R1/L + mu11/S1 ) / DUMMY
	 X11(4)=-( L*R1 + X11(2) ) * ( N*( 1.D0/L + mu11)/R1   ) / DUMMY
         DUMMY = DSQRT( X11(1)**2 + X11(2)**2 + X11(3)**2  + X11(4)**2)
	 X11(1)=X11(1)/DUMMY
	 X11(2)=X11(2)/DUMMY
	 X11(3)=X11(3)/DUMMY
	 X11(4)=X11(4)/DUMMY

	 DUMMY = ( (1.D0-S1)/L ) * ( (1.D0+A)*R1/L + mu12/s1 ) - (N/R1)*( 1.D0/L + mu12 )*( R1/L + mu12/S1 )
	 X12(1)=0.D0
	 X12(2)=1.D0
	 X12(3)=- ( (1.D0+A)*R1/L + mu12/S1 ) / DUMMY
	 X12(4)=- ( N*( 1.D0/L + mu12 )/R1  ) /DUMMY
         DUMMY = DSQRT( X12(1)**2 + X12(2)**2 + X12(3)**2 + X12(4)**2)
	 X12(1)=X12(1)/DUMMY
	 X12(2)=X12(2)/DUMMY
	 X12(3)=X12(3)/DUMMY
	 X12(4)=X12(4)/DUMMY


	 X14(1)=0.D0
	 X14(2)=0.D0
	 X14(3)= ( R1/L + mu14/S1 ) / ( (1.D0-S1)/L )
	 X14(4)=1.D0
         DUMMY = DSQRT( X14(1)**2 + X14(2)**2 + X14(3)**2 + X14(4)**2)
	 X14(1)=X14(1)/DUMMY
	 X14(2)=X14(2)/DUMMY
	 X14(3)=X14(3)/DUMMY
	 X14(4)=X14(4)/DUMMY

	 X13(1)=0.D0
	 X13(2)=0.D0
	 X13(3)=1.D0
	 X13(4)=N* ((1.D0-S1)/L) / ( N*R1/L + N*mu13/S1)
         DUMMY = DSQRT( X13(1)**2 + X13(2)**2 + X13(3)**2 + X13(4)**2)
	 X13(1)=X13(1)/DUMMY
	 X13(2)=X13(2)/DUMMY
	 X13(3)=X13(3)/DUMMY
	 X13(4)=X13(4)/DUMMY

	 X04(1)=0.D0
	 X04(2)=0.D0
	 X04(3)= ( R0/L + mu04/S0 ) / ( (1.D0-S0)/L )
	 X04(4)=1.D0
         DUMMY = DSQRT( X04(1)**2 + X04(2)**2 + X04(3)**2 + X04(4)**2)
	 X04(1)=X04(1)/DUMMY
	 X04(2)=X04(2)/DUMMY
	 X04(3)=X04(3)/DUMMY
	 X04(4)=X04(4)/DUMMY



! unit length constraints of the coefficients
	FB(1)= (PAR(9)**2 + PAR(10)**2 + PAR(11)**2) -1.0d0
	FB(2)= (PAR(12)**2 + PAR(13)**2 + PAR(14)**2) -1.d0

! boundary values
	FB(3) = U0(1) - ( 0.0D0 + EPS0*(  PAR( 9)*X01(1) + PAR(10)*X02(1) + PAR(11)*X03(1)   ) )
	FB(4) = U0(2) - ( 0.0D0 + EPS0*(  PAR( 9)*X01(2) + PAR(10)*X02(2) + PAR(11)*X03(2)   ) )
	FB(5) = U0(3) - ( R0    + EPS0*(  PAR( 9)*X01(3) + PAR(10)*X02(3) + PAR(11)*X03(3)   ) )
	FB(6) = U0(4) - ( S0    + EPS0*(  PAR( 9)*X01(4) + PAR(10)*X02(4) + PAR(11)*X03(4)   ) )


	FB(7) = U1(1) - ( 0.0D0 + EPS1*(  PAR(12)*X11(1) + PAR(13)*X12(1) + PAR(14)*X14(1)   ) )
	FB(8) = U1(2) - ( 1.0D0 + EPS1*(  PAR(12)*X11(2) + PAR(13)*X12(2) + PAR(14)*X14(2)   ) )
	FB(9) = U1(3) - ( R1    + EPS1*(  PAR(12)*X11(3) + PAR(13)*X12(3) + PAR(14)*X14(3)   ) )
	FB(10)= U1(4) - ( S1    + EPS1*(  PAR(12)*X11(4) + PAR(13)*X12(4) + PAR(14)*X14(4)   ) )

! use explicit formulas for eigenvalues
!	FB(24)=PAR(28)-mu01
!	FB(25)=PAR(29)-mu02
!	FB(26)=PAR(30)-mu03
!	FB(27)=PAR(31)-mu11
!	FB(28)=PAR(32)-mu12




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




      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS





