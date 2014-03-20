      subroutine rs(rx_max,ry_max, imax, jmax, x_temp, y_temp, rx, ry)
	  implicit none
	  integer imax, jmax, i, j, k
	  double precision rx(imax-1, jmax-1), rx_max, ry_max
	  double precision ry(imax-1, jmax-1)
	  double precision mx
	  double precision a, b, c
	  double precision xe, xw,xs,xn
	  double precision xen, xes, xwn, xws
	  double precision x(imax+1, jmax+1), y(imax+1, jmax+1)
	  double precision x_temp(imax+1, jmax+1), y_temp(imax+1, jmax+1)
	  double precision x_sigma, y_sigma, x_yita, y_yita
C..  check the residual and store it into r
      do i=2, imax
	      do j=2, jmax
              ! the derivatives
			  x_sigma=(x_temp(i+1,j)-x_temp(i-1,j))/2
			  y_sigma=(y_temp(i+1,j)-y_temp(i-1,j))/2
			  x_yita=(x_temp(i,j+1)-x_temp(i,j-1))/2
			  y_yita=(y_temp(i,j+1)-y_temp(i,j-1))/2
			  ! the coefficient: a, b, c
			  a=x_yita**2+y_yita**2
			  b=x_sigma*x_yita+y_sigma*y_yita
			  c=x_sigma**2+y_sigma**2
			  ! calculate the new value for x and update x_temp
			  xe=x_temp(i+1,j)
			  xw=x_temp(i-1,j)
			  xn=x_temp(i,j+1)
			  xs=x_temp(i,j-1)
			  xen=x_temp(i+1,j+1)
			  xes=x_temp(i+1,j-1)
			  xwn=x_temp(i-1,j+1)
			  xws=x_temp(i-1,j-1)				  
			  rx(i-1,j-1)=abs(2*(a+c)*x_temp(i,j)-a*(xe+xw)+
     >                   (0.5*b*(xen-xes+xws-xwn)-c*(xn+xs)))
			  ! calculate the new value for y and update y_temp
			  xe=y_temp(i+1,j)
			  xw=y_temp(i-1,j)
			  xn=y_temp(i,j+1)
			  xs=y_temp(i,j-1)
			  xen=y_temp(i+1,j+1)
			  xes=y_temp(i+1,j-1)
			  xwn=y_temp(i-1,j+1)
			  xws=y_temp(i-1,j-1)				  
			  ry(i-1,j-1)=abs(2*(a+c)*y_temp(i,j)-a*(xe+xw)+
     >                   (0.5*b*(xen-xes+xws-xwn)-c*(xn+xs)))	 
		  end do
	  end do
C..  choose the biggest value of residual and give it to r_max
      ! x
      mx=rx(1,1)
      do i=1, imax-1
	      do j=1, jmax-1
		      if (rx(i,j)>=mx) then
			      mx=rx(i,j)
			  else 
			      continue
			  end if
		  end do
	  end do
	  rx_max=mx
	  ! y
      mx=ry(1,1)
      do i=1, imax-1
	      do j=1, jmax-1
		      if (ry(i,j)>=mx) then
			      mx=ry(i,j)
			  else 
			      continue
			  end if
		  end do
	  end do
	  ry_max=mx
	  return
	  end 
		      