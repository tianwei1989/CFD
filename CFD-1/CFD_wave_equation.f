      Program solve_wave_equation
	  implicit none
C.. define variable
      integer i,j, k, n,imax, count_no
	  integer case_no
	  real dx, lx(100), dt, t
	  real pi
	  real cfl
	  real c
	  real u(100,100)
	  real u_ana(100)
	  character(len=100) filename

C.. parameter initialization for control
	  !parameter (case_no=2)
	  !parameter (cfl=0.3)
	  !parameter (n=1)
	  parameter (imax=41)
C.. parameter initialization
      parameter(dx=1.0)
	  parameter (c=1)
	  parameter (pi=3.1415)
C.. variable initialization
      i=1
	  j=1
	  lx(1)=0.0
	  t=0.0
	  count_no=1
C.. read the parameters
      write (6,*) "Case_no=?"
      read (5,*) case_no
      write (6,*) "n=?"
      read (5,*) n
      write (6,*) "cfl=?"
      read (5,*) cfl
C.. This is the 8th assignment of CFD course
      write (6,*) 'This is the 8th assignment of CFD course'
C.. space discretization
      do i=2,41
	      lx(i)=lx(i-1)+dx
	  end do
C.. give initial value
      do i=1,41
	      u(i,1)=sin(2*n*pi*lx(i)/40)
	  end do
C.. use CFL and dx to calculate dt
      dt=cfl*dx/c  
C.. Routine for different schemes
      select case (case_no)
C.. 1--Lax Wendroff one step scheme
	      case (1)
		      write (6,*) 'Entrance for the Lax-Wendroff one step method'
              call lax_wendroff (u,lx,imax,c,dx,dt,t,j)
C.. 2--MacCormack Two Step scheme		  
		  case (2)
		      write (6,*) 'Entrance for the MacCormack Two Step scheme'
              call Mac_Cormack (u,lx,imax,c,dx,dt,t,j)	  
C.. 3--Second order upwind scheme		  
		  case (3)
		      write (6,*) 'Entrance for the Second order upwind scheme'  
	          call second_order_upwind (u,lx,imax,c,dx,dt,t,j)	   
C.. 4--4th order Runge-Kutta scheme		  
		  case (4)
		      write (6,*) 'Entrance for 4th order Runge-Kutta scheme'
              call runger_kutta (u,lx,imax,c,dx,dt,t,j)			  
C.. default to output error information	  
		  case default
		      write (6,'(A)') 'can not recognize the selector i'  
      end select
C.. delete me
C.. Output data 
      write (filename,*) "count",".txt"
      open (1, file='results.txt', status='unknown')
	  !open (1, file='wave_equ.txt', status='unknown', access='append')
	  if (case_no==1) then
	      write(1,*) 'results from Lax-Wendroff one step scheme'
	  else if (case_no==2) then
		  write(1,*) 'results from MacCormack two step scheme'
	  else if (case_no==3) then
		  write(1,*) 'results from second order upwind scheme'
	  else if (case_no==4) then
		  write(1,*) 'results from 4th order Runge-Kutta scheme'
	  else
          write (1,*) 'Error with outputting data'
      end if
	  write (1,*) 'cfl=', cfl, 'n=',n, 'time step size=', dt
	  write (1,2001) 
2001  format (5x,'X',10x,'Numerical solution',5x,'Analytical solution')
	  do i=1,41
	      u_ana(i)=sin(2*n*pi*(lx(i)-18)/40)
	      write (1,2002) lx(i),u(i,j ), u_ana(i)
2002  format (f13.10,6x, f13.10, 10x, f13.10)
	  end do
	  count_no=count_no+1
      End program solve_wave_equation