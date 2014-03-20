      subroutine GSL_SOR(imax,jmax,x,y,x_temp,y_temp)
	  implicit none
C..  use Gauss_Seidel line iteration
C..  use the new value immediately when available
C..  line sweeps vertically from west to east
      integer imax, jmax, i, j, k
	  double precision a, b, c, sor 
	  double precision xe, xw,xs,xn
	  double precision xen, xes, xwn, xws
	  double precision x(imax+1, jmax+1), y(imax+1, jmax+1)
	  double precision x_temp(imax+1, jmax+1), y_temp(imax+1, jmax+1)
	  double precision x_temp1(imax+1, jmax+1)
	  double precision y_temp1(imax+1, jmax+1)
	  double precision x_sigma, y_sigma, x_yita, y_yita
C.. define the variable for TDMA
      double precision ax(imax-1), bx(imax-1), cx(imax-1), dx(imax-1)
	  double precision ay(imax-1), by(imax-1), cy(imax-1), dy(imax-1)
C..  define a local matrix to store the residual, r_o is to output residual
      double precision rx(imax-1, jmax-1), ry(imax-1, jmax-1)
	  double precision rx_max, ry_max
	  double precision rx_o(100000), ry_o(100000)
C..  define SOR-- relaxing coefficient. Default 2.
      write (6, *) 'please input SOR coefficient, default 2'
	  write (6, *) 'maximal 1.5'
	  read (5, *) sor
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
			      ! store the value of x_temp and v_temp at 
				  ! x_temp1 and y_temp1
				  x_temp1(i,j)=x_temp(i,j)
				  y_temp1(i,j)=y_temp(i,j)
			      ! the derivatives
			      x_sigma=(x_temp(i+1,j)-x_temp(i-1,j))/2
			      y_sigma=(y_temp(i+1,j)-y_temp(i-1,j))/2
                  x_yita=(x_temp(i,j+1)-x_temp(i,j-1))/2
			      y_yita=(y_temp(i,j+1)-y_temp(i,j-1))/2
			      a=x_yita**2+y_yita**2
			      b=x_sigma*x_yita+y_sigma*y_yita
			      c=x_sigma**2+y_sigma**2
			      ! calculate the new variable of x to simplify formula
			      xe=x_temp(i+1,j)
				  xw=x_temp(i-1,j)
				  xn=x_temp(i,j+1)
				  xs=x_temp(i,j-1)
				  xen=x_temp(i+1,j+1)
				  xes=x_temp(i+1,j-1)
				  xwn=x_temp(i-1,j+1)
				  xws=x_temp(i-1,j-1)
				  ! calculate ax bx cx for X
                  ax(j-1)=c
				  bx(j-1)=-2*(a+c)
				  cx(j-1)=c
				  ! calculate dx for X
				  dx(j-1)=0.5*b*(xen-xes+xws-xwn)-a*(xe+xw)
				  if (j==2) then
				      dx(j-1)=dx(j-1)-c*x_temp(i,j-1)				  
				  else if (j==jmax) then
				      dx(j-1)=dx(j-1)-c*x_temp(i,j+1)	
				  else
				      continue
				  end if				  
				  ! calculate the new variable of x to simplify formula
			      xe=y_temp(i+1,j)
				  xw=y_temp(i-1,j)
				  xn=y_temp(i,j+1)
				  xs=y_temp(i,j-1)
				  xen=y_temp(i+1,j+1)
				  xes=y_temp(i+1,j-1)
				  xwn=y_temp(i-1,j+1)
				  xws=y_temp(i-1,j-1)
				  ! calculate ay by cy for Y
                  ay(j-1)=c
				  by(j-1)=-2*(a+c)
				  cy(j-1)=c
				  ! calculate dy for Y
				  dy(j-1)=0.5*b*(xen-xes+xws-xwn)-a*(xe+xw)
				  if (j==2) then
				      dy(j-1)=dy(j-1)-c*y_temp(i,j-1)				  
				  else if (j==jmax) then
				      dy(j-1)=dy(j-1)-c*y_temp(i,j+1)	
				  else
				      continue
				  end if
			  end do
		      call sy(jmax-1,ax, bx, cx, dx)
			  call sy(jmax-1,ay, by, cy, dy)
              do j=2,jmax
                  x_temp(i,j)=dx(j-1)
                  y_temp(i,j)=dy(j-1)
              end do
			  ! update the x_temp and y_temp with SOR
			  do j=2, jmax
			      x_temp(i,j)=x_temp1(i,j)+sor*(x_temp(i,j)-x_temp1(i,j))
				  y_temp(i,j)=y_temp1(i,j)+sor*(y_temp(i,j)-y_temp1(i,j))
              end do		  
		  end do
          call rs(rx_max,ry_max, imax, jmax, x_temp, y_temp, rx, ry)
		  k=k+1
		  rx_o(k)=rx_max
	      ry_o(k)=ry_max          
          write (6, *)"residual, No.", k , rx_max, ry_max
      end do
C..  give all the value of x_temp and y_temp to x and y for output
      do i=2, imax
	      do j=2, jmax
		      x(i,j)=x_temp(i,j)
			  y(i,j)=y_temp(i,j)
		  end do
	  end do
C..  output residual
	  open (2, file='gsl_sor_rs.dat', status='unknown')
	  write (2, *) 'title=residual'
	  do i=1,k
	      write (2, *) i, rx_o(i), ry_o(i)
	  end do
	  write (6, *) "k is", k
	  return
	  end	  
		  