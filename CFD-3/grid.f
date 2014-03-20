      program grid
C..   This codes is to solve the 11th CFD homework
C..   Both uniform and stretched grids are developed
      implicit none
C..   Define Variables
C..   mesh size in x and y
      integer imax, jmax, case_no, i,j,k
C..   matrix used to store the coordinates
      double precision x(100,100), y(100,100), x_temp, y_temp 
C..   coefficient in stretched grids formulation
      double precision a 
	  double precision pi
C..   geometry 
      double precision sita, beta, l_1, l_2, l_3, x_t
      double precision y_t, x_c, y_c, r_c, h_i
C..   initialize the data
      parameter (imax=51, jmax=25)
	  parameter (pi=3.1415926)
	  parameter (sita=-22.33/360*2*pi,beta=1.21/360*2*pi)
	  parameter (l_1=4.74, l_2=5.84, l_3=11.56, x_t=5.84)
      parameter (y_t=1.37, x_c=5.78, y_c=4.11, r_c=2.74, h_i=3.52)
C..  input the case no
      write (6, *) 'please input the case no'
	  write (6, *) '1: uniform; 2: stretched'
	  read (5, *) case_no
C..	  calculate the coordinates
      select case (case_no)
	      ! uniform grids
	      case (1)
		      ! coordinates in x
		      do i=1, imax+1
			      do j=1, jmax+1
			          x(i,j)=(i-1)*l_3/imax
				  end do
			  end do
			  ! calculate y and the coordinates as x 
              do i=1, imax+1
			      x_temp=x(i,1)
			      if (x_temp>=0 .and. x_temp<=l_1) then
				      do j=1, jmax+1
				          y_temp=tan(sita)*x_temp+h_i
						  y(i,j)=(j-1)*y_temp/jmax
					  end do
				  else if (x_temp>=l_1 .and. x_temp<=l_2) then
				      do j=1, jmax+1
				          y_temp=-1*sqrt(r_c**2-(x_temp-x_c)**2)+y_c
						  y(i,j)=(j-1)*y_temp/jmax
					  end do
				  else 
				      do j=1, jmax+1
				          y_temp=tan(beta)*(x_temp-x_t)+y_t
						  y(i,j)=(j-1)*y_temp/jmax
					  end do
                  end if
              end do				  

		  ! stretched grids
		  case (2)
		      ! input the coefficient for equation generating stretched grids
		      write (6,*) 'input equation coefficient a'
		      read (5,*) a
		      ! coordinates in x
		      do i=1, imax+1
			      do j=1, jmax+1
			          x(i,j)=(i-1)*l_3/imax
				  end do
			  end do
			  ! calculate y and the coordinates as x 
              do i=1, imax+1
			      x_temp=x(i,1)
			      if (x_temp>=0 .and. x_temp<=l_1) then
				      y_temp=tan(sita)*x_temp+h_i
				      do j=1, jmax+1
						  y(i,j)=y_temp-(real(j-1)/jmax)**a*y_temp
					  end do
				  else if (x_temp>=l_1 .and. x_temp<=l_2) then
				      y_temp=-1*sqrt(r_c**2-(x_temp-x_c)**2)+y_c
					  do j=1, jmax+1
						  y(i,j)=y_temp-(real(j-1)/jmax)**a*y_temp
					  end do
				  else 
				      y_temp=tan(beta)*(x_temp-x_t)+y_t
					  do j=1, jmax+1
						  y(i,j)=y_temp-(real(j-1)/jmax)**a*y_temp				  
					  end do
                  end if
              end do	  
		  case default
		     write (6, *) 'Error with the case input'
          end select
C..   output data in Tecplot format
      open (1,file='results.plt', status='unknown')
C..   write the title of the figure
	  write (1,*) "variable= x, y"
	  write (1,*) "ZONE I=52, J=26, F=POINT"
C..   output data line by line
	  do j=1, jmax+1
	      do i=1, imax+1
	          write (1,*) x(i,j), y(i,j)
		  end do
	  end do
C..  end the program
	  end 	  