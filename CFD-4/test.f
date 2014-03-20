      subroutine reader(imax, jmax, case_no, a, r)
	  implicit none
	  integer imax, jmax
	  character case_no
	  double precision a, r
	  open (1,file='input.cfd', status='unknown')
	  read (1,*)
	  read (1,*) imax
	  read (1,*)
	  read (1,*) jmax
	  read (1,*)
	  read (1,*) case_no
	  read (1,*)
	  read (1,*) a
	  read (1,*)
	  read (1,*) r
	  
	  write (6, *) imax, jmax, case_no, a
	  end 