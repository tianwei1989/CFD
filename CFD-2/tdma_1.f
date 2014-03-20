      subroutine sy(iu,bb, dd, aa, cc)
	  implicit none
	  integer il, iu, lp, i, j
	  double precision bb(iu), dd(iu), aa(iu), cc(iu), r
	  il=1
	  lp=il+1
	  do i=lp, iu
	      r=bb(i)/dd(i-1)
		  dd(i)=dd(i)-r*aa(i-1)
		  cc(i)=cc(i)-r*cc(i-1)
	  end do
	  
	  cc(iu)=cc(iu)/dd(iu)
	  do i=lp,iu
	      j=iu-i+il
		  cc(j)=(cc(j)-aa(j)*cc(j+1))/dd(j)
	  end do
	  return
	  end