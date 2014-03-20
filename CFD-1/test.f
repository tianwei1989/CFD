      program CFD_course_wave_equation

c..   CFD Course, MEN 614

c..   This is to calculate the wave equation in different schemes. 
c..   Zha 02/02/2002
c..   The wave equation is: du/dt + du/dx =0
c..   The initial solution is u(x,0)=sin(2*n*pi*x/40.), 0<=x<=40.
c..   The analytical solution will be u(x,t) = sin(2*n*pi*(x-t)/40.)
c..   Periodic condition will be used at the two ends x=0 and x=40.
c..   Chose 41 grid point mesh with Dx=1 and compute to t=18.
c..   Solve the problem with n=1,3, and CFL=1.0, 0.6,0.3

      implicit none

      integer n                 ! frequency in initial solution
      integer ngrid             ! number of grid points
      double precision x0,xl    ! x-coordinates of starting and ending points
      double precision cfl      ! CFL number
      integer iteration         ! number of iteration for upper limit
      double precision t_stop   ! time to stop the calculation
      double precision time     ! marching time level

      parameter(ngrid=41,x0=0.,xl=40.,iteration=100, t_stop=18.)

      double precision u_n(ngrid)        ! Solution at time level n
      double precision u_nplus1(ngrid)   ! Solution at time level n+1
      double precision u_analytical      ! analytical solution

      double precision x                 ! x-coordinates
      double precision dx                ! grid spacing delta x
      double precision dt                ! time step
      double precision c                 ! wave speed of the linear wave eq.
      double precision pi                ! 3.14159

      integer i,j            ! loop index
      integer scheme         ! scheme name
                             ! =1, Lax scheme
                             ! =2, Lax-Wendroff Scheme



      write(6,*)'Please input n and CFL number'
      read(5,*)n,cfl
      write(6,*)'What scheme to use?'
      write(6,*)'If to use Lax scheme, please input 1'
      read(5,*)scheme

c..   For the linear wave equation in this problem, c=1

      c=1.0
      pi = 3.14159

c..   give the initial solution

      dx = (xl-x0)/float(ngrid-1)

      do i = 1,ngrid
         x    = x0 + float(i-1)*dx  
         u_n(i) = sin(2.*float(n)*pi*x/40.)
      end do

c..   determine time step

      dt = CFL*dx/c

c..   Start Time Marching

      IF (scheme.eq.1) then

c..      Lax Scheme
         
         time =0.0
         
         DO i = 1, iteration

            do j = 2, ngrid-1

            u_nplus1(j) = 0.5*(u_n(j+1) + u_n(j-1))
     >                   - c*dt/dx*0.5*(u_n(j+1) - u_n(j-1))

            end do

c..         implement periodic boundary condition for starting point

            u_nplus1(1) = 0.5*(u_n(2) + u_n(ngrid-1))
     >                   - c*dt/dx*0.5*(u_n(2) - u_n(ngrid-1))

c..         implement periodic boundary condition for ending  point

            u_nplus1(ngrid) = u_nplus1(1)
            

            time = time + dt

c..      save u_nplus1 to u_n for next iteration
            
            do j = 1, ngrid
               u_n(j) = u_nplus1(j) 
            end do

            if(time.ge.t_stop) go to 111

         END DO

      ELSE IF (scheme.eq.2) then

c..      Other scheme will be implemented

         continue

      END IF

c...  Output results

 111  continue

c..   Output on screen

      if(scheme.eq.1)write(6,*)'Lax Scheme Results at time=', time
      write(6,*)'n=',n
      write(6,*)'CFL=',cfl

      write(6,113)
 113  format(13x,'x',20x,'u_cfd',20x,'u_analytical')

c..   Output to a file

      open(1,file='waveeq_results')

      if(scheme.eq.1)write(1,*)'Lax Scheme Results at time=', time
      write(1,*)'n=',n
      write(1,*)'CFL=',cfl

      write(1,113)

      do i = 1, ngrid
         x    = x0 + float(i-1)*dx  
         u_analytical = sin(2.*float(n)*pi*(x-time)/40.)
         write(6,*)x,u_nplus1(i),u_analytical
         write(1,*)x,u_nplus1(i),u_analytical
      end do

      STOP
      END
