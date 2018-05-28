!***************************************************************
!*          Sorting an array with the Shell method             *
!* ----------------------------------------------------------- *
!* REFERENCE:                                                  *
!*      "NUMERICAL RECIPES by W.H. Press, B.P. Flannery,       *
!*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
!*       University Press, 1986".                              *
!* ----------------------------------------------------------- *
!***************************************************************
subroutine SORT2(A, MAX_VALUE, N)

integer :: N, MAX_VALUE 
real    :: A(MAX_VALUE)     !Table to be sorted

! call sorting subroutine
  call SHELL(N,A)

  return

end

!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!*            by the Shell-Mezgar method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*            N          size of table ARR           *
!*          ARR          table to be sorted          *
!* OUTPUT:                                           *
!*            ARR   table sorted in ascending order  *
!*                                                   *
!* NOTE: The Shell method is a N^3/2 routine and can *
!*       be used for relatively large arrays.        * 
!*****************************************************         
subroutine SHELL(N,ARR)
real, parameter :: ALN2I=1./0.69314718, TINY=1.E-5
  integer:: n, i, j, k, l, m, nn, lognb2
  real ARR(N), t
  LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
  m=n
  do nn=1,LOGNB2
    m=m/2; k=n-m
    do j=1,k
      i=j
10    continue
      l=i+m
      if(ARR(l).LT.ARR(i)) then
        t=ARR(i)
        ARR(i)=ARR(l)
        ARR(l)=t
        i=i-m
        if(i.GE.1) GOTO 10
      end if
    end do
  end do
  return
end

!write table of size N to standard output
subroutine TWRIT(N,ARR)
  integer:: n, i
  real ARR(N)

  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
10 FORMAT(10F6.1)
end

!end of file sort2.f90
