      program plate
C.. This program is used to calculate the velocity distribution within two plates,
C.. of which one is suddenly moving at u.
      implicit none
C.. define the variable
      double precision u(100), u0(100), u_temp(100,100)
	  double precision dx, l,lx(100), t, dt, tmax, nu
	  double precision ub, u_init
	  double precision cfl
	  integer imax, i, j, k, case_no
C..  initialize the value
      write (6,*) 'please input the case no.'
	  write (6,*) '1 for EXPLICIT; 2 for IMPLICIT'
	  read (5, *) case_no
	  write (6,*) 'please input the time step size'
	  write (6,*) 'dt=0.00200'
	  read (5,*) dt
      i=0
	  j=0
	  k=0
	  t=0.00000
	  l=0.0400
	  dx=0.00100
	  tmax=1.08000
	  ub=40.00000
	  u_init=40.00000
	  nu=0.000217
C.. Discrete the space
      imax=nint(l/dx)
	  do i=1,imax+1
	      lx(i)=(i-1)*dx
	  end do
C.. give initial value
      u0(1)=ub
      do i=2,imax+1
	      u0(i)=0
	  end do   
C.. initialize the t
	  t=dt
C..  set boundary conditions
      u0(1)=ub
      u0(imax+1)=0
	  u(1)=ub
      u(imax+1)=0	  
C.. main routine for computation
      select case (case_no)
	      case (1)
		      write (6,*) "Upwind explicit schemes"
              call explict_upwind (u,u0,imax,dx,dt,t,tmax,nu,ub,u_temp)
          case (2)
              write (6,*) "Upwind implicit schemes"
              call implict_upwind (u,u0,imax,dx,dt,t,tmax,nu,ub,u_temp)
          case default
              write (6,*) "Error with selecting the schemes"
	      end select
C.. output data
	  open (1, file='results.txt', status='unknown')
	  write (6, *) 'simulation time is',tmax,'seconds'
	  write (6, *),'time step size is',dt, 'seconds'
	  write (6, *) 'computation length is',l, 'm'
	  write (1, 2001)
2001  format (10x,'NO.', 10x, 'X', 23x, 'Velocity') 
      do i=1,imax+1
	      write (1, *) i,lx(i),u0(i)
	  end do
	  read (5, *)
	  end
			  