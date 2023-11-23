!==============================================================================!
  pure subroutine Gaussian_Elimination(Math, n, a, b, x, invertible)
!------------------------------------------------------------------------------!
!>  Gaussian elimination with scaled row pivoting.
!------------------------------------------------------------------------------!
!   From: Numerical Analysis:                                                  !
!         The Mathematics of Scientific Computing                              !
!         D.R. Kincaid & E.W. Cheney                                           !
!         Brooks/Cole Publ., 1990                                              !
!         Section 4.3                                                          !
!         https://web.ma.utexas.edu/CNA/cheney-kincaid/f90code/CHP07/gauss.f90 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Math_Type), intent(in)    :: Math        !! parent class
  integer,          intent(in)    :: n           !! number of unknows
  real,             intent(inout) :: a(n,n)      !! system matrix
  real,             intent(inout) :: b(n)        !! right hand side vector
  real,             intent(inout) :: x(n)        !! solution vectoru
  logical,          intent(inout) :: invertible  !! true if matrix invertible
!-----------------------------------[Locals]-----------------------------------!
  real,    allocatable :: s(:)
  integer, allocatable :: p(:)
  real                 :: r, rmax, smax, sum, z
  integer              :: i, k, j, pk
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Math)
!==============================================================================!

  allocate(s(n))
  allocate(p(n))

  invertible = .true.

  do i = 1, n
    p(i) = i
    smax = 0.0
    do j = 1, n
      smax = max(smax, abs(a(i,j)))
    end do
    s(i) = smax
  end do

  do k = 1, n-1
    rmax = 0.0
    do i = k, n
      r = abs(a(p(i),k)) / s(p(i))
      if (r .gt. rmax) then
        j = i
        rmax = r
      endif
    end do

    pk   = p(j)
    p(j) = p(k)
    p(k) = pk

    do i = k+1, n
      z = a(p(i),k) / a(p(k),k)
      a(p(i),k) = z
      do j = k + 1, n
        a(p(i),j) = a(p(i),j) - z * a(p(k),j)
      end do
      if(maxval(abs(a(p(i),k+1:n))) < PICO) then
        invertible = .false.
        return
      end if
    end do
  end do


  do k = 1, n-1
    do i = k + 1, n
      b(p(i)) = b(p(i)) - a(p(i),k) * b(p(k))
    end do
  end do

  do i = n, 1, -1
    sum = b(p(i))
    do j = i + 1, n
      sum = sum - a(p(i),j) * x(j)
    end do
    x(i) = sum / a(p(i),i)
  end do

  deallocate(s)
  deallocate(p)

  end subroutine
