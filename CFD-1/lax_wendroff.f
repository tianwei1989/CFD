      subroutine lax_wendroff (u,lx,imax,c,dx,dt,t,j)
C..  this is the subroutine for lax wendroff method
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
C.. main entrance for lax wendroff method	  
	  do while (t<=18)
		  do i=1,41
		 ! Special treatment of point 1
		      if (i==1) then
				  u(i,j+1) = u(i,j)
     >                       - c*dt/2/dx*(u(i+1,j)-u(imax-1,j)) 				     
     > 					     + c*c*dt*dt/2/dx/dx*(u(i+1,j)-2*u(i,j) 
     >                       + u(imax-1,j))   
          ! Special treatment of point 41	 
			  else if (i==41) then
				  u(i,j+1) = u(i,j)-c*dt/2/dx*(u(2,j)-u(i-1,j))      
     >					     + c*c*dt*dt/2/dx/dx*(u(i-1,j)-2*u(i,j)    
     >				         + u(2,j))   
	      ! the regular point
			  else 
				  u(i,j+1) = u(i,j)-c*dt/2/dx*(u(i+1,j)-u(i-1,j))     
     >					    + c*c*dt*dt/2/dx/dx*(u(i+1,j)-2*u(i,j)     
     >					    + u(i-1,j))
			  end if
          end do					  
		  j=j+1
		  t=t+dt
	  end do
      return
      end 