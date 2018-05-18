program lab1
	integer, parameter :: N = 8
	real :: X(N)=(/0.0, 0.2, 0.5, 0.7, 1.0, 1.3, 1.7, 2.0/)
	real :: Y(N)=(/1.000, 1.1487, 1.4142, 1.6245, 2.0000, 2.4623, 3.2490, 4.0000/) 
	real :: B(N), C(N), D(N)
	call findRootUsingSpline(N, X, Y)
end

subroutine findRootUsingSpline(N, X, Y)
	integer :: N, I
	real :: X(N), Y(N), B(N), C(N), D(N)
	real :: X1, X2, S
	real :: splineRes, iterTemp, delta, tmp
	X1 = X(1)
	X2 = X(N)
	S = 1

	call SPLINE(N, X, Y, B, C, D) 
	do while (S > 0.00001)
		delta = (X2 - X1) * 0.5
		tmp = X1 + delta
		S = SEVAL(N, tmp, X, Y, B, C, D) + equation(tmp)
		write(*,"(f10.4, f10.4)") tmp, S
		if (S > 0) then
			X2 = tmp
		else
			X1 = tmp
			S = -S 
		end if
	end do
end subroutine
	
real function equation(x)
	real x
	equation = 5 * x - 3
	return
end

SUBROUTINE SPLINE (N,X,Y,B,C,D)
!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------
      IMPLICIT REAL *4 (A-H,O-Z)
      DIMENSION B(N),C(N),D(N),X(N),Y(N)
      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 50

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO 10 I = 2,NM1
      D(I) = X(I+1)-X(I)
      B(I) = 2.D0*(D(I-1)+D(I))
      C(I+1) = (Y(I+1)-Y(I))/D(I)
      C(I) = C(I+1)-C(I)
   10 CONTINUE

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0
      IF (N.EQ.3) GO TO 15
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

   15 DO 20 I = 2,N
      T = D(I-1)/B(I-1)
      B(I) = B(I)-T*D(I-1)
      C(I) = C(I)-T*C(I-1)
   20 CONTINUE

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO 30 L = 1,NM1
      I = N-L
      C(I) = (C(I)-D(I)*C(I+1))/B(I)
   30 CONTINUE

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
      DO 40 I = 1,NM1
      B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
      D(I) = (C(I+1)-C(I))/D(I)
      C(I) = 3.D0*C(I)
   40 CONTINUE
      C(N) = 3.D0*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.D0
      D(1) = 0.D0
      B(2) = B(1)
      C(2) = 0.D0
      D(2) = 0.D0
      RETURN
      END

 FUNCTION SEVAL (N,U,X,Y,B,C,D)
!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     SEVAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------
      REAL *4 B(N),C(N),D(N),X(N),Y(N),U,DX,SEVAL
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I)) GO TO 10
      IF (U.LE.X(I+1)) GO TO 30
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 20

!     SPLINE EVALUATION

   30 DX = U-X(I)
      SEVAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      RETURN
      END