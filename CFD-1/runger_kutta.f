      subroutine runger_kutta (u,lx,imax,c,dx,dt,t,j)
C..  this is the subroutine for Runger-Kutta method	
      implicit none
	  real u(100,100)
	  real lx(100)
	  integer imax
	  real c,dx,dt
	  real u_temp1(100), u_temp2(100),u_temp3(100)
	  real R1(100), R2(100), R3(100), R4(100)
	  integer i,j,k
	  real t
C..  variable initialization
      i=1
	  j=1
	  k=1
	  t=dt
C.. main entrance for Runger-Kutta method
      do while (t<=18)
      ! first step
	      ! calculate R1
		  do i=1,imax
		      if (i==1) then
			      R1(i)=-1*c*(u(i+1,j)-u(imax-1,j))/(2*dx)
			  else if (i==41) then 
			      R1(i)=-1*c*(u(2,j)-u(i-1,j))/(2*dx)
			  else
			      R1(i)=-1*c*(u(i+1,j)-u(i-1,j))/(2*dx)
			  end if 
		  end do
	      ! calculate u_temp1
          do i=1,imax
			  u_temp1(i)=u(i,j)+0.5*dt*R1(i)
          end do 
	  ! Second step
          ! calculate R2
		  do i=1,imax
		      if (i==1) then
			      R2(i)=-1*c*(u_temp1(i+1)-u_temp1(imax-1))/(2*dx)
			  else if (i==41) then 
			      R2(i)=-1*c*(u_temp1(2)-u_temp1(i-1))/(2*dx)
			  else
			      R2(i)=-1*c*(u_temp1(i+1)-u_temp1(i-1))/(2*dx)
			  end if 
		  end do   
          ! Calculate u_temp2
          do i=1,imax
			  u_temp2(i)=u(i,j)+0.5*dt*R2(i)
          end do 	  
	  ! Third step
	      ! calculate R3
		  do i=1,imax
		      if (i==1) then
			      R3(i)=-1*c*(u_temp2(i+1)-u_temp2(imax-1))/(2*dx)
			  else if (i==41) then 
			      R3(i)=-1*c*(u_temp2(2)-u_temp2(i-1))/(2*dx)
			  else
			      R3(i)=-1*c*(u_temp2(i+1)-u_temp2(i-1))/(2*dx)
			  end if 
		  end do 
          ! calculate u_temp3
          do i=1,imax
			  u_temp3(i)=u(i,j)+dt*R3(i)
          end do 
 	      ! calculate R4
		  do i=1,imax
		      if (i==1) then
			      R4(i)=-1*c*(u_temp3(i+1)-u_temp3(imax-1))/(2*dx)
			  else if (i==41) then 
			      R4(i)=-1*c*(u_temp3(2)-u_temp3(i-1))/(2*dx)
			  else
			      R4(i)=-1*c*(u_temp3(i+1)-u_temp3(i-1))/(2*dx)
			  end if 
		  end do
      ! Forth Step
	      do i=1,imax
		      u(i,j+1)=u(i,j)+dt/6*(R1(i)+2*R2(i)+2*R3(i)+R4(i))
		  end do
      t=t+dt
	  j=j+1
      end do	
      return
      end 	  