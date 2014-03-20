      subroutine set_bnd(u,ub,imax)
	  implicit none
	  double precision u(100)
	  real ub
	  integer imax
C.. Set the boundary conditions
	  u(imax+1)=0
	  u(1)=ub
	  return
	  end
	  