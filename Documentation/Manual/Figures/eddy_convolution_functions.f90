program functions

implicit none

integer, parameter :: n =  50
real, parameter    :: lx    = 1.0
real, parameter    :: sigma = 0.5
real, parameter    :: pi    = 3.14159265359

integer :: i
real    :: x

do i = 0, n

  x = -lx / 2.0 + i * lx / n

  print *, x,                                                      &
           x,                                                      &
           1.0 / (sigma *sqrt(PI+PI))*exp(-0.5*((x-0)/sigma)**2),  &
           x   / (sigma *sqrt(PI+PI))*exp(-0.5*((x-0)/sigma)**2)

end do

end program

