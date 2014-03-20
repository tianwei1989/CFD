      program burger_equation
	  implicit none
	  integer imax, ileft, iright, tmax
	  character (100) filename
	  parameter (imax=51, ileft=11, iright=40)
	  parameter (tmax=50)
	  double precision u_ini_l, u_ini_r
	  parameter (u_ini_l=1.00000)
	  parameter (u_ini_r=0.00000)
	  double precision u_o(imax, tmax)
	  double precision u_temp(imax),u (imax)
      integer i, j, k, case_no, io
      double precision cfl, dt,dx, lx, t
	  double precision fp, fe, fw, ue, uw, up, ae, aw
	  double precision up_bar, uw_bar, fp_bar, fw_bar
	  double precision r_n(imax), r_1(imax), r_2(imax), r_3(imax)
	  double precision u1(imax), u2(imax), u3(imax)
C..  input cfl number
      write (6, *) 'input cfl number'
	  read (5, *) cfl
C..  determine dx and dx
      lx=1.000
	  dx=lx/(imax-1)
	  dt=cfl*dx
	  t=0.00000
C..  initialize u_o
      do i=1, imax
	      do j=1, tmax
		     u_o(1,j)=0.00000
		  end do
	  end do
C..  initialization
      do i=1, imax
          if (i<=ileft) then
              u(i)=	u_ini_l
			  u_temp(i)= u_ini_l
			  u_o(i,1)=	u_ini_l
		  else 
		      u(i)=	u_ini_r
			  u_temp(i)= u_ini_r
			  u_o(i,1)=	u_ini_r
		  end if
	  end do
C..  choose scheme
      write (6, *) 'choose scheme'
	  write (6, *) '1--Lax first order scheme'
	  write (6, *) '2--Lax-Wendroff One Step Scheme'
	  write (6, *) '3--Lax-MacCormack Two Step Scheme'
	  write (6, *) '4--4th Runge-Kutta Scheme'
	  read (5, *) case_no
C..  calculate routine
      select case (case_no)
	      case (1)
		      k=0
		      write (6, *) '1--Lax first order scheme'
			  do while (k<=tmax)
				  do i=2, imax-1
					  ue=u_temp(i+1)
					  uw=u_temp(i-1)			  
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  u(i)=0.5*(ue+uw)-dt/dx*0.5*(fe-fw)
				  end do
				  ! update u_temp
				  do i=2, imax-1
					  u_temp(i)=u(i)
				  end do
				  !time and k forward one step
				  t=t+dt
				  k=k+1
				  !give value to u_o for output
				  do i=1, imax-1
					  u_o(i,k+1)=u(i)
				  end do
			  end do		  
		  case (2)
		      write (6, *) '2--Lax-Wendroff One Step Scheme'
			  k=0
			  do while (k<=tmax)
                  do i=2, imax-1
					  ue=u_temp(i+1)
					  uw=u_temp(i-1)
					  up=u_temp(i)
                      fp=0.5*up**2					  
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  ae=0.5*(up+ue)
					  aw=0.5*(up+uw)
					  u(i)=up-dt/dx*0.5*(fe-fw)+0.5*(dt/dx)**2
     >                      *(ae*(fe-fp)-aw*(fp-fw))
	              end do
				  ! update u_temp
				  do i=2, imax-1
					  u_temp(i)=u(i)
				  end do
				  !time and k forward one step
				  t=t+dt
				  k=k+1
				  !give value to u_o for output
				  do i=1, imax-1
					  u_o(i,k+1)=u(i)
				  end do
			  end do		  
		  case (3)
		      write (6, *) '3--Lax-MacCormack Two Step Scheme'
			  k=0
			  do while (k<=tmax)
                  do i=2, imax-1
					  ue=u_temp(i+1)
					  uw=u_temp(i-1)
					  up=u_temp(i)
                      fp=0.5*up**2					  
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  ae=0.5*(up+ue)
					  aw=0.5*(up+uw)
					  up_bar=up-dt/dx*(fe-fp)
					  uw_bar=uw-dt/dx*(fp-fw)
					  fp_bar=0.5*up_bar**2
					  fw_bar=0.5*uw_bar**2
					  u(i)=0.5*(up+up_bar-dt/dx*(fp_bar-fw_bar))
				  end do
				  ! update u_temp
				  do i=2, imax-1
					  u_temp(i)=u(i)
				  end do
				  !time and k forward one step
				  t=t+dt
				  k=k+1
				  !give value to u_o for output
				  do i=1, imax-1
					  u_o(i,k+1)=u(i)
				  end do
			  end do					  
		  case (4)
		      write (6, *) '4--4th Runge-Kutta Scheme'
			  k=0
			  do while (k<=tmax)
			      do i=2, imax-1
					  ue=u_temp(i+1)
					  uw=u_temp(i-1)
					  up=u_temp(i)
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  r_n(i)=-(fe-fw)/2/dx	
                      u1(i)=up+dt/2*r_n(i)	
                  end do
                  u1(1)=u_ini_l
                  u1(imax)=u_ini_r
				  ! second step
			      do i=2, imax-1
					  ue=u1(i+1)
					  uw=u1(i-1)
					  up=u_temp(i)
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  r_1(i)=-(fe-fw)/2/dx	
                      u2(i)=up+dt/2*r_1(i)	
                  end do
                  u2(1)=u_ini_l
                  u2(imax)=u_ini_r
				  ! Third step
			      do i=2, imax-1
					  ue=u2(i+1)
					  uw=u2(i-1)
					  up=u_temp(i)
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  r_2(i)=-(fe-fw)/2/dx	
                      u3(i)=up+dt*r_2(i)	
                  end do
                  u3(1)=u_ini_l
                  u3(imax)=u_ini_r
				  ! Forth step
			      do i=2, imax-1
					  ue=u3(i+1)
					  uw=u3(i-1)
					  up=u_temp(i)
					  fe=0.5*ue**2
					  fw=0.5*uw**2
					  r_3(i)=-(fe-fw)/2/dx	
                      u(i)=up+dt/6*(r_n(i)+2*r_1(i)+2*r_2(i)+r_3(i))	
                  end do
				  ! update u_temp
				  do i=2, imax-1
					  u_temp(i)=u(i)
				  end do
				  !time and k forward one step
				  t=t+dt
				  k=k+1
				  !give value to u_o for output
				  do i=1, imax-1
					  u_o(i,k+1)=u(i)
				  end do
			  end do					  
      end select
C..  output data
      write (filename, '(A6, I1.1, A4)') 'result', case_no, '.dat'
	  filename=trim(filename)
      open (1, file=filename,status='unknown')
	  write (1, '(A7, 1X, 10A5)') 'Title=', 'Lax'
	  do i=1, imax
		  write (1, *) i,(u_o(i,j), j=1,tmax)
      end do		    
C..  end program
      end
	  
	  
 			  