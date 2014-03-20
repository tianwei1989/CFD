      subroutine explict_upwind (u,u0,imax,dx,dt,t,tmax,nu,ub,u_temp)
	  implicit none
      double precision u(100)
	  double precision u0(100), u_temp(100,100)
	  integer imax, i, j
	  double precision dx, dt,t,tmax,nu,ub
	  double precision small
	  j=2
	  do while (t<=tmax)
		  do i=2,imax
			  u(i)=u0(i)+nu*dt/dx/dx*(u0(i+1)-2*u0(i)+u0(i-1))
		  end do
		  do i=2,imax
			  u0(i)=u(i)
		  end do 
          t=t+dt		  
	  end do 
	  return
	  end