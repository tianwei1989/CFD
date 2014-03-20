      subroutine GSP(imax,jmax,x,y,x_temp,y_temp)
	  implicit none
C..  use Gauss_Seidel point iteration
C..  use the new value immediately when available
      integer imax, jmax, i, j, k
	  double precision a, b, c 
	  double precision xe, xw,xs,xn
	  double precision xen, xes, xwn, xws
	  double precision x(imax+1, jmax+1), y(imax+1, jmax+1)
	  double precision x_temp(imax+1, jmax+1), y_temp(imax+1, jmax+1)
	  double precision x_sigma, y_sigma, x_yita, y_yita
C..  define a local matrix to store the residual, r_o is to output residual
      double precision rx(imax-1, jmax-1), ry(imax-1, jmax-1)
	  double precision rx_max, ry_max
	  double precision rx_o(100000), ry_o(100000)
c.. check the residual and store the initial value to r_o
      call rs(rx_max,ry_max, imax, jmax, x_temp, y_temp, rx, ry)
	  rx_o(1)=rx_max
	  ry_o(1)=ry_max
	  write (6, *)"residual for x is " , rx(2,2)
C..  set j=2 to continue store residual
      k=1
C..  main routine for calculation
      do while (rx_max>1e-12 .or. ry_max>1e-12)
	  !do while (k<=10000)
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
				  x(i,j)=0.5*a/(a+c)*(xe+xw)-0.5/(a+c)
     >                   *(0.5*b*(xen-xes+xws-xwn)-c*(xn+xs))
			      x_temp(i,j)=x(i,j)
			     ! calculate the new value for y and update y_temp
			      xe=y_temp(i+1,j)
				  xw=y_temp(i-1,j)
				  xn=y_temp(i,j+1)
				  xs=y_temp(i,j-1)
				  xen=y_temp(i+1,j+1)
				  xes=y_temp(i+1,j-1)
				  xwn=y_temp(i-1,j+1)
				  xws=y_temp(i-1,j-1)				  
				  y(i,j)=0.5*a/(a+c)*(xe+xw)-0.5/(a+c)
     >                   *(0.5*b*(xen-xes+xws-xwn)-c*(xn+xs))
			      y_temp(i,j)=y(i,j)			  			  		  	 
			  end do
		  end do
		  ! check the residual
		  call rs(rx_max,ry_max, imax, jmax, x_temp, y_temp, rx, ry)
		  k=k+1
		  rx_o(k)=rx_max
	      ry_o(k)=ry_max          
          write (6, *)"residual, No.", k , rx_max, ry_max	  
	  end do
	  open (2, file='gsp_rs.dat', status='unknown')
	  write (2, *) 'title=residual'
	  do i=1,k
	      write (2, *) i, rx_o(i), ry_o(i)
	  end do
	  write (6, *) "k is", k
	  return
	  end
      