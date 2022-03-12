!==============================================================================!
  program Stretch
!------------------------------------------------------------------------------!
!   Places the nodes on the line defined with local block position             !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: n, i, case
  real              :: x0, delx, t, dt, ddt, pr, xi
  real              :: x(0:1000)
  real              :: w
  character(len=20) :: file1 
!==============================================================================!

  x0=0.0
 
  print *, 'Distance: '
  read(*,*)  delx 
  print *, 'Weigth: '
  read(*,*)  w
  print *, 'Number of cells: '
  read(*,*)  n

  x(0) = x0

  !-------------------------!
  !   Linear distribution   !
  !-------------------------!
  if(w  > 0.0) then
    ddt = ( 2.0*(1.0-w) ) / ( 1.0*n*(1.0*n-1.0)*(1.0+w) )
    t=0.0
    do i=1,n
      dt=1.0/(1.0*n)+(1.0*i-0.5*(1.0*n+1)) * ddt
      t=t+dt
      x(i) = x0 + t*delx
    end do

  !-----------------------------!
  !   Hyperbolic distribution   !
  !-----------------------------!
  else
    case = 0
    if     ((w  >  -0.5).and.(w <=  -0.25)) then
      pr = 1.0 - abs(0.5 - abs(w))    
      case = 1
    else if((w >=  -0.75).and.(w  <  -0.5)) then
      pr = 1.0 - abs(0.5 - abs(w))            
      case = 2
    else
      pr = -w
      case = 3 
    endif

    do i=1,n
      if(case .eq. 1) xi = -1.0*(1.0*i)/(1.0*n)
      if(case .eq. 2) xi =  1.0 - 1.0*(1.0*i)/(1.0*n)
      if(case .eq. 3) xi = -1.0 + 2.0*(1.0*i)/(1.0*n)
      if    (case .eq. 1) then
        x(i) = x0 - (tanh(xi*atanh(pr))/pr)*delx
      elseif(case .eq. 2) then
        x(i) = x0 + delx - (tanh(xi*atanh(pr))/pr)*delx
      elseif(case .eq. 3) then
        x(i) = x0 + 0.5*(1.0+tanh(xi*atanh(pr))/pr)*delx
      endif 
    end do

  endif

  !--------------------------------!
  !   Print out the distribution   !
  !--------------------------------!
  print *, ' number:   node:                  cell:'
  do i=n,1,-1
    write(*,'(I6,F12.6,A10)') i+1, x(i),              ' +-------+'
    write(*,'(I6,A23,F12.6)') i,  '|   o   | ', 0.5*(x(i)+x(i-1))
  end do
  write(*,'(I6,F12.6,A10)') i+1, x(0),              ' +-------+'
  print *, 'wall:    //////////////////////////////'

  file1 = 'fileout'

  open(3,file=file1) 
  
  do i = 1, n
  
   write(3,*) 0.5*(x(i)+x(i-1)) 

  end do

  close (3)
 
  end program 

