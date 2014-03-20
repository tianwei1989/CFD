      subroutine Mac_Cormack (u,lx,imax,c,dx,dt,t,j)
C..  this is the subroutine for the Mac Cormack method
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
C.. main entrance for Mac Cormack method
      do while (t<=18)
      ! first step-prediction
	      do i=1,imax
              if (i==imax) then 
                  u_temp(i,j+1)=u(i,j)-c*dt/dx*(u(2,j)-u(i,j))
              else 				  
			      u_temp(i,j+1)=u(i,j)-c*dt/dx*(u(i+1,j)-u(i,j))
              end if 
          end do
	  ! second step prediction	
          do i=1,imax
              if (i==1) then
			      u(i,j+1)=0.5*(u(i,j)+u_temp(i,j+1)
     >				    -c*dt/dx*(u_temp(i,j+1)-u_temp(imax-1,j+1)))	
              else
                  u(i,j+1)=0.5*(u(i,j)+u_temp(i,j+1)
     >				    -c*dt/dx*(u_temp(i,j+1)-u_temp(i-1,j+1)))
              end if 	 
	      end do
		  t=t+dt
		  j=j+1
      end do
      return 
      end