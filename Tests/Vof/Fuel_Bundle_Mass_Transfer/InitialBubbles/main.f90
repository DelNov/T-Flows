  PROGRAM MAIN
!
!      0 <= X,Y <= 0.066
!  -0.06 <= Z <= 0.121
!
  IMPLICIT NONE
  INTEGER::n, i, seedsize
  REAL(8):: rand
  REAL(8):: x, y, z, r
  REAL(8):: x0, x1, alx, y0, y1, aly, z0, z1, alz
  REAL(8):: r0, r1, alr
  integer,allocatable::seed(:)

  x0 = 0.0d0;   x1=1.280d-2; alx = x1-x0
  y0 = 0.0d0;   y1=1.280d-2; aly = y1-y0
  z0 = 0.0d0+2d-3; z1=3.936d-02 - 2d-3; alz = z1-z0  ! don't put near inlet and outlet
  r0 = 0.2d-3;  r1=0.6d-3;  alr = r1-r0

  n = 500  !1000

  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)

  write(*,*)"seedsize =",seedsize
  write(*,*)"seed =",seed
  
  OPEN(10, file='ellipsoid_parameters.ini')
    WRITE(10,*)n
    do i=1,n
      call random_number(rand)
      x = x0 + rand * alx
      call random_number(rand)
      y = y0 + rand * aly
      call random_number(rand)
      z = z0 + rand * alz
      call random_number(rand)
      r = r0 + rand * alr
      WRITE(10,*)
      write(10,'(3E13.5)')r, r, r
      write(10,'(3E13.5)')x, y, z
    end do
  CLOSE(10)

  STOP
  END
