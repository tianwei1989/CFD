      subroutine bc(u,v,p,p_pre,imax,jmax,dx,dy)
	  implicit none
C..   **************************************************
C..                       Note
C..   1. Four corners are neglected
C..   2. recently only pressure B.C. is set here
C..   Author: Wei, Tian
C..   Date: 3/20/2014
C..   **************************************************
	  integer imax, jmax
	  double precision dx, dy
	  double precision u(imax+2, jmax+2), v(imax+2, jmax+2)
	  double precision p(imax+2, jmax+2), p_re(imax+1, jmax+1)
      ! set the last one of p_re to be zero
	      p_re(imax+1, jmax+1)=0
	  ! update the pressure of the fluid cell
	  do i=2, imax+1
	      do j=2, jmax+1
		      p(i,j)=p(i,j)+p_re(i,j)
		  end do
	  end do
	  ! update the pressure of the solid cell
	  ! south and north
	  do i=2, imax+1
	      j=1
		  p(i,j)=p(i,j+1)
		  j=jmax+2
		  p(i,j)=p(i,j-1)
	  end do
	  ! west and east
	  do j=2, jmax+1
	      i=1
		  p(i,j)=p(i+1,j)
	      i=imax+2
		  p(i,j)=p(i-1,j)
	  end do
	  
	  return
	  end 
		  
	      
	  
	  
	  