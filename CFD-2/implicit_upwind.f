      subroutine implict_upwind (u,u0,imax,dx,dt,t,tmax, nu,ub,u_temp)
	  implicit none
	  integer imax, i,j
	  double precision u(100)
	  double precision u0(100), u_temp(100,100)
	  double precision a(imax-1), b(imax-1), c(imax-1), d(imax-1)
	  double precision x(imax-1)
	  double precision dx, dt,t,tmax,nu,ub, psi
	  integer n
	  n=imax-1
C.. calculate psi
      psi=nu*dt/dx/dx
	  write (6,*) "coefficient for equation is", psi
C.. Start do iteration
	  do while (t<=tmax)
C..  upper diagonal coefficient
		  do i=1,n
			  a(i)=psi
		  end do
C..  main diagonal coefficient
		  do i=1,n
			  b(i)=-1*(2*psi+1)
		  end do
C..  lower diagonal coefficient
		  do i=1,n
			  c(i)=psi
		  end do

C.. d constant vector
		  do i=1,n
			  if (i==n) then
				  d(i)=-u0(imax)-psi*u0(imax+1)
			  else if (i==1) then
				  d(i)=-u0(i+1)-psi*u0(i)
			  else 
				  d(i)=-u0(i+1)
			  end if
		  end do
		  !call tdma(n,a,b,c,d,x)
		  call sy(n,a,b, c, d)
C.. update the u0
		  do i=2, imax
			  u0(i)=d(i-1)
		  end do 
C.. iteration continue
		  t=t+dt
	  end do
	  return
	  end 