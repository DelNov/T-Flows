!==============================================================================!
  subroutine Brent(Fun, x0, tol, res, low_bnd, upp_bnd,  &
                   detail, lhs_lin, lhs_fun, rhs)
!------------------------------------------------------------------------------!
! Find root of a function using Brent's method.                                !
!                                                                              !
! Inputs:                                                                      !
!                                                                              !
! Fun:                     subroutine to evaluate (jump_condition)             !
! x0:                      initial guess                                       !
! tol:                     error tolerance                                     !
! low_bnd, upp_bnd:  bounds for the root                                       !
! detail:                  output result of iteration if set to 1              !
! lhs_lin, lhs_funs, rhs:  arguments for jump condition                        !
!------------------------------------------------------------------------------!
! Description of algorithm:                                                    !
!                                                                              !
! Firstly, this subroutine tries to find x1 and x2 such that sign of f(x1)     !
! and f(x2) are not equal.                                                     !
!                                                                              !
! Secondly, it finds a root of function f(x) given intial bracketing interval  !
! [a,b] where f(a) and f(b) must have opposite signs. At a typical step we     !
! have three points a, b, and c such that f(b)f(c)<0, and a may coincide with  !
! c. The points a, b, and c change during the algorithm, and the root always   !
! lies in either [b,c] or [c, b]. The value b is the best approximation to     !
! the root and a is the previous value of b.                                   !
!                                                                              !
! The iteration uses following selection of algorithms:                        !
!                                                                              !
! * when bracket shrinks reasonablly fast,                                     !
!   - linear interporation if a == b                                           !
!   - quadratic interporation if a != b and the point is in the bracket;       !
!                                                                              !
! * othrwise:                                                                  !
!   - bisection.                                                               !
!                                                                              !
! Based on https://github.com/yoki/Optimization                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    subroutine  Fun(x, fval0, lhs_lin, lhs_fun, rhs)
    implicit none
    real :: x
    real:: fval0
    real :: lhs_lin, lhs_fun, rhs
    end subroutine Fun
  end interface
!---------------------------------[Arguments]----------------------------------!
  real,    intent(in) :: x0
  real,    intent(in) :: tol
  real                :: res
  real,    intent(in) :: low_bnd
  real,    intent(in) :: upp_bnd
  integer, intent(in) :: detail
  real                :: lhs_lin
  real                :: lhs_fun
  real                :: rhs
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, iter, exitflag, disp
  real    :: a , b , olda, oldb, fa, fb
  real    :: dx  ! change in bracket
  real    :: sgn
  real    :: c, diff,e, fc, p, q, r, s, tol1,xm,tmp 
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: D       = selected_real_kind(p=15, r=307) ! precision
  integer, parameter :: MAXITER = 40
  integer, parameter :: IMAX    = 100  ! maximum number of iteration
  real,    parameter :: SQRT2   = sqrt(2.0)! change in dx
  real,    parameter :: EPS     = epsilon(a)
  ! values from Numerical Recipe
!==============================================================================!

  a  = x0  ! lower bracket
  b =  x0 ! upper bracket
  olda = a  
  oldb = b 
  exitflag = 0  ! flag to see we found the bracket
  call Fun(x0, sgn, lhs_lin, lhs_fun, rhs) ! sign of initial guess

  disp = detail  

  ! Set initial change dx
  if (abs(x0)<0.00000002) then 
    dx = 1.0/50.0
  else
    dx = 1.0/50.0 * x0
  end if

  if (disp == 1) then
    print *, 'Search for initial guess for Brents method'
    print *, 'find two points whose sign for f(x) is different '
    print *, 'x1 searches downwards, x2 searches ' //  &
             'upwards with increasing increment'
    print *, ' '
    print *, '   i         x1         x2         f(x1)          f(x2)'
  end if


  ! Main loop to extend a and b
  do iter = 1, MAXITER
    call Fun(a,fa,lhs_lin, lhs_fun, rhs)
    call Fun(b,fb,lhs_lin, lhs_fun, rhs)

    if (disp == 1) write(*,"(1I4,4F14.7)") iter, a, b, fa, fb

    ! check if sign of functions changed or not
    if ( (sgn >= 0 ) .and.  (fa <= 0) ) then  ! sign of a changed 
      ! use a and olda as bracket
      b = olda
      exitflag = 1
      exit
    else if  ( (sgn <= 0 ) .and.  (fa >= 0  ) ) then ! sign of b changed
      b = olda
      exitflag = 1
      exit
    else if  ( (sgn >= 0 ) .and.  (fb <= 0  ) ) then ! sign of a changed
      a = oldb
      exitflag = 1
      exit
    else if  ( (sgn <= 0 ) .and.  (fb >= 0  ) ) then ! sign of a changed
      a = oldb
      exitflag = 1
      exit
    end if

    ! Update boundary
    olda = a
    oldb = b
    a  = a - dx
    b  = b+ dx
    dx = dx * SQRT2

    ! Boundary check
    if (a < low_bnd ) a = low_bnd + tol
    if (b > upp_bnd ) b = upp_bnd - tol

  end do

  if (exitflag /=  1 ) then 
    print *, ' Error (brent2) : Proper initial value ' //  &
             'for Brents method could not be found'
    print *, ' Change initial guess and try again. '
    print *, ' You might want to try disp = 1 option too'
    print *, 'i              x1            x2          fx1              fx2'
    write(*,"(1i4,4F12.7)") iter, a, b, fa, fb
    print *, '  press any key to abort the program'
    read(*,*) 
    stop
  else if (disp == 1) then
    print *, '  Initial guess was found.'
    print *, ''
  end if

  exitflag = 0

  !----------------------!
  !   Intialize values   !
  !----------------------!
  c = b
  call Fun(a, fa, lhs_lin, lhs_fun, rhs)
  call Fun(b, fb, lhs_lin, lhs_fun, rhs)
  fc = fb

  ! Check sign
  if ( (fa>0. .and. fb>0. )  .or.  (fa>0. .and. fb>0. )) then
    print *,  'Error (brent.f90): Root must be bracked by two imputs'
    write(*, "('f(x1) = ', 1F8.4, '   f(x2) = ', 1F8.4)") fa,fb
    print *, 'press any key to halt the program'
    read(*,*)
    stop
  end if

  if (disp == 1 ) then 
    print *, 'Brents method to find a root of f(x)'
    print *, ' '
    print *, '  i           x          bracketsize            f(x)'
  end if

  !--------------------!
  !   Main iteration   !
  !--------------------!
  do i = 1, IMAX

    ! Rename c and adjust bounding interval if both a(=b) and c are same sign
    if ((fb > 0.  .and. fc > 0) .or. (fb <0. .and. fc < 0. ) ) then 
      c = a
      fc = fa
      e = b-a
      diff = e
    end if

    ! If c is better guess than b, use it.
    if (abs(fc) < abs(fb) ) then 
      a  = b
      b  = c
      c  = a
      fa = fb
      fb = fc
      fc = fa
    end if

    ! Convergence check
    tol1=2.0* EPS * abs(b) + 0.5*tol
    xm = 0.5 * (c - b)
    if (abs(xm) < tol1 .or. fb == 0.0 )  then
      exitflag = 1
      exit
    end if

    if (disp == 1) then 
      tmp = c-b
      write(*,"('  ', 1I2, 3F16.6)") i, b, abs(b-c), fb
    end if

    ! Try inverse quadratic interpolation
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then 
      s = fb/fa
        if (abs(a - c) < EPS) then 
        p = 2.0 *xm * s
        q = 1.0  - s
      else
        q = fa/fc
        r = fb/fc
        p = s * (2.0 * xm * q * (q -r ) - (b - a) * (r - 1.0))
        q = (q - 1.0 ) * (r - 1.0) * (s - 1.0) 
      end if

      ! Accept if q is not too small to stay in bound
      if (p > 0.0) q = -q
      p = abs(p)
      if (2.0 * p < min(3.0 * xm * q - abs(tol1* q), abs(e *q))) then 
        e = D
        diff = p / q
      else   ! interpolation failed. use bisection
        diff= xm 
        e = D
      end if
    else  ! quadratic interpolation bounds moves too slowly, use bisection
      diff = xm
      e = D
    end if

    ! Update last bound
    a  = b
    fa = fb

    ! Move the best guess
    if (abs(D) > tol1) then 
      b = b + diff
    else
      b = b + sign(tol1, xm)
    end if

    ! Evaluate new trial root
    call Fun(b, fb, lhs_lin, lhs_fun, rhs)
  end do

  !------------------------------!
  !   Case for non convergence   !
  !------------------------------!
  if (exitflag /= 1 ) then
    print *, 'Error (brent.f90) :  convergence was not attained'
    print *, 'Initial value:'
    write(*,"(4f10.5)" )   a, b, fa, fb
    print *, ' '
    print *, 'final value:'
    write(*,"('x = '  ,1F6.4, ':                f(x1) = ' ,  1F6.4  )" )   &
         b,  fb
  else if( disp == 1) then
    print *, 'Brents method was converged.'
    print *, ''
  end if

  res = b

  return

end subroutine
