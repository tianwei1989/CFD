      subroutine second_order_upwind (u,lx,imax,c,dx,dt,t,j)	
C..  this is the subroutine for the  second order upwind method
      implicit none
	  real u(100,100)
	  real lx(100)
	  integer imax
	  real c,dx,dt
	  real u_temp(100,100)
	  integer i,j,k
	  real t
C..  variable initialization
      i=1
	  j=1
	  k=1
	  t=dt
C.. main entrance for second order upwind method
      do while (t<=18)
	      !first step-predictor
          do i=1,imax
		      if (i==1) then
			      u_temp(i,j+1)=u(i,j)-c*dt/dx*(u(i,j)-u(imax-1,j))
			  else 
			      u_temp(i,j+1)=u(i,j)-c*dt/dx*(u(i,j)-u(i-1,j))
			  end if 
		  end do
		  !second step-corrector
		  do i=1,imax
		      if (i==1) then 
			      u(i,j+1)=0.5*(u(i,j)+u_temp(i,j+1)-c*dt/dx*(u_temp(i,j+1)
     >                     -u_temp(imax-1,j+1))-c*dt/dx*(u(i,j)
     >                     -2*u(imax-1,j)+u(imax-2,j)))
			  else if (i==2) then 
			      u(i,j+1)=0.5*(u(i,j)+u_temp(i,j+1)-c*dt/dx*(u_temp(i,j+1)
     >                     -u_temp(i-1,j+1))-c*dt/dx*(u(i,j)
     >                     -2*u(i-1,j)+u(imax-1,j)))			  
			  else
			      u(i,j+1)=0.5*(u(i,j)+u_temp(i,j+1)-c*dt/dx*(u_temp(i,j+1)
     >                     -u_temp(i-1,j+1))-c*dt/dx*(u(i,j)
     >                     -2*u(i-1,j)+u(i-2,j)))			 
			  end if
		  end do
		  t=t+dt
		  j=j+1
	  end do
	  return
	  end