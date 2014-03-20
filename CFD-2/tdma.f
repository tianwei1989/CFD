      subroutine tdma(n,a,b,c,d,x)
	  implicit none
      integer n
      double precision a(n), c(n)
      double precision b(n), d(n)
	  double precision x(n)
	  !  --- Local variables ---
	  integer i
	  double precision q
      !  --- Elimination ---
	  write (6,*) "n in the iteration is ", n
      do i = 2,n
         q = a(i)/b(i - 1)
	     write (6,*) "q in the iteration is ", q
         b(i) = b(i) - c(i - 1)*q
         d(i) = d(i) - d(i - 1)*q
		 write (6,*) "dn is", d(i)
      end do
	  !write (6,*) "bn is ", b(n), "dn is", d(n)
      ! --- Backsubstitution ---
      q = d(n)/b(n)
	  !write (6,*) "q in the iteration is ", q
      x(n) = q
      do i = n - 1,1,-1
         q = (d(i) - c(i)*q)/b(i)
         x(i) = q
      end do
      return
      end