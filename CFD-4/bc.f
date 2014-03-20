      subroutine set_bc(sig,imax,jmax,r,x,y,x_temp,y_temp)
	  implicit none
C..  this subroutine is used to set boundary conditions for x and x_temp
C..  the inner and outer circle is separated evenly
C..  define the variables
      integer imax, jmax, i, j, k, sig
	  double precision r, a
	  double precision x(imax+1, jmax+1), y(imax+1, jmax+1)
	  double precision x_temp(imax+1, jmax+1), y_temp(imax+1, jmax+1)
C..   calculate angle	  
	  double precision alpha
	  double precision pi
      parameter (pi=3.1415926)
C..  set boundary conditions for outer circle, north, which is (i,jmax+1)
      do i=1, imax+1
	      j=jmax+1
		  x(i,j)=cos(-1*real(i-1)/imax*(2*pi))*(5+1)*r
		  y(i,j)=sin(-1*real(i-1)/imax*(2*pi))*(5+1)*r
		  x_temp(i,j)=cos(-1*real(i-1)/imax*(2*pi))*(5+1)*r
		  y_temp(i,j)=sin(-1*real(i-1)/imax*(2*pi))*(5+1)*r
	  end do
C..  set boundary conditions for inner circle, south, which is (i,1)
      do i=1, imax+1
	      j=1
		  x(i,j)=cos(-1*real(i-1)/imax*(2*pi))*r
		  y(i,j)=sin(-1*real(i-1)/imax*(2*pi))*r
		  x_temp(i,j)=cos(-1*real(i-1)/imax*(2*pi))*r
		  y_temp(i,j)=sin(-1*real(i-1)/imax*(2*pi))*r
	  end do
C..  set boundary for west and east, which is (1,j) 
      if (sig==1) then  
		  do j=1, jmax+1
			  ! west
			  i=1
			  x(i,j)=r+real(j-1)/jmax*5*r
			  y(i,j)=0
			  x_temp(i,j)=r+real(j-1)/jmax*5*r
			  y_temp(i,j)=0
			  ! east
			  i=imax+1
			  x(i,j)=r+real(j-1)/jmax*5*r
			  y(i,j)=0	
			  x_temp(i,j)=r+real(j-1)/jmax*5*r
			  y_temp(i,j)=0			  
		  end do
	  else if (sig==2) then
	      write (6, *) 'input value of a'
		  write (6, *) 'larger than 2 required'
		  read (5, *) a
	      do j=1, jmax+1
		      !west
			  i=1
			  x(i,j)=r+(real(j-1)/jmax)**a*5*r
			  y(i,j)=0
			  x_temp(i,j)=r+(real(j-1)/jmax)**a*5*r
			  y_temp(i,j)=0
		      ! east
			  i=imax+1
			  x(i,j)=r+(real(j-1)/jmax)**a*5*r
			  y(i,j)=0
			  x_temp(i,j)=r+(real(j-1)/jmax)**a*5*r
			  y_temp(i,j)=0
		  end do
	  else
	      write (6, *) 'Error!!!'
	  end if
	  return
	  end