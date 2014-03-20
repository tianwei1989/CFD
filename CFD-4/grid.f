      program elliptic_grid
	  implicit none
C..  this program aims to finish the assignment of Project 1 in CFD homework
C..  The differential equation methods, or elliptic scheme , is used.
C..  Outer layer is at 5 times of radius of the inner circle
C..  Gauss-Seidel point, Gauss-Seidel line sweep and Gauss-Seidel sweep with SOR are used
C..  Both uniform and stretched grids can be generated by this program
C..  bc.f specifies BC
C..  gsp.f is for Gauss-Seidel Point method
C..  gsl.f is for Gauss-Seidel line sweep method
C..  gsp_sor.f is for Gauss-Seidel Line sweep with SOR
C..  rs.f is for checking the residual of each iteration
C..  tdma.f is the Tridiagonal Matrix Algorithm which is used by gsl.f
C..         and gsl_sor.f
C..  Author: Wei Tian
C..  C#11345638
C..  Date:3/05/2014

C..  define the case selector for scheme
      character(100) case_no
C..  define imax and jmax fraction in x and y, which ends up in imax+1 and jmax+1 points
      integer imax, jmax, i, j, k, sig
C..  define the value of imax and jmax
	  parameter (imax=30, jmax=30)
C..  define the radius:r
      double precision r
      parameter (r=1.00000)
C..  define the variable to store x and y coordinates and corresponding temporal one
      double precision x(imax+1, jmax+1), y(imax+1, jmax+1)
	  double precision x_temp(imax+1, jmax+1), y_temp(imax+1, jmax+1)
C..  initialize the variable
      do i=1, imax+1
	      do j=1, jmax+1
		      x(i,j)=0
			  y(i,j)=0
			  x_temp(i,j)=0
			  y_temp(i,j)=0
		  end do
	  end do
C.. input uniform method or stretched method
      write (6, *) 'please input uniform or stretched grids'
	  write (6, *) '1:uniform; 2: stretched'
	  read (5, *) sig
C..  input the scheme want to use
      write (6, *) 'input the scheme'
	  write (6, *) 'm-1 : Gauss Seidel Point Iteration'
	  write (6, *) 'm-2 : Gauss Seidel Line Iteration'
	  write (6, *) 'm-3 : Gauss Seidel Line Iteration with SOR'
	  read (5, *) case_no
C..  call the subroutine to give the boundary conditions
      call set_bc(sig,imax,jmax,r,x,y,x_temp,y_temp)
C.. main routine for calculating
      	  select case (case_no)
		      case ('m-1')
			      call GSP(imax,jmax,x,y,x_temp,y_temp)
			  case ('m-2')
			      call GSL(imax,jmax,x,y,x_temp,y_temp)
			  case ('m-3')
			      call GSL_SOR(imax,jmax,x,y,x_temp,y_temp)
			  case default
			      write (6, *) 'Error!!!!!!'
		  end select
C..  output the result in tecplot format
      open (1, file='result.plt', status='unknown')
      write (1,*) "variable= x, y"
	  write (1,*) "ZONE I=",imax+1,',','J=',jmax+1,'F=POINT' 
      do j=1, jmax+1
          do i=1, imax+1
              write (1, *) x(i,j), y(i,j)
          end do
      end do
      write (6, *)"X(2,2) " , x_temp(2,11)	
      read (5,*)	  
	  end 