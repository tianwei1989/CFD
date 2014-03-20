      program cfd
C..   this program calculates steady-state Navier-Stokes equation
C..   2-D case in X and Y
C..   lid-driven cavity whose upper jar is moving at speed u0
C..   Finite Volume Method is used in discretizing the space
C..   SIMPLE scheme is used for solving pressure correction
C..   Mass balance equation is additionally used to determine the iteration loop
C..   Gauss-Seidel is used to do calculation
C..   Results are output in Tecplot format
C..   Author: Wei Tian
C..   Date: 3/14/2014
      implicit none
C..   **********************************************************************************
C..   **********************************************************************************
C..   define the variable
      double precision lx, ly, dx,dy, nu
	  !* imax and jmax mean how many inner points
	  integer imax, jmax, i, j, k
	  parameter (i=10, j=10, lx=1.00000, ly=1.00000)
	  !* the dimension should be imax+2 and jmax+2 
	  double precision u(imax+2, jmax+2), v(imax+2, jmax+2)
	  double precision p(imax+2, jmax+2), p_re(imax+1, jmax+1)
	  !* physical properties
	  parameter (nu=1.5e-5)
	  !***********************************
	  ! the velocity running at of the jar
	  double precision u0
	  parameter (u0=1.00000)
	  !***********************************
C..   *************************************************************************************
C..   *************************************************************************************
C..   mesh generation. fixme: later on, algebraic and elliptic method can be used
      dx=lx/imax
	  dy=ly/jmax
C..   initialization for u
      do i=1, imax+2
	      do j=1, jmax+2
		      u(i,j)=0
		  end do
	  end do
	  do i=1, imax+2
	      j=jmax+2
		  u(i, j)=u0
	  end do
C..   initialization for v, p, 
      do i=1, imax+1
	      do j=1, jmax+1
              v(i, j)=0.00000
              p(i, j)=0.00000
          end do
      end do
C..   initialization for p_re	  
      do i=1, imax
	      do j=1, jmax
              p_re(i, j)=0.00000
          end do
      end do

	  end do
      