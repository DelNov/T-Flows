!==============================================================================!
  subroutine Cg(Nat, Acon, Aval, x, b, miter, tol)
!------------------------------------------------------------------------------!
!   Note: This is an alternative algorithm and I am honestly not sure where    !
!         I have found it any more, but it avoids one "if" block during the    !
!         iterations.  It gives exactly the same result as the original Cg.    !
!         This verision was implemented in SFS package, and I wanted to make   !
!         sure it is the same as the original Cg.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),    target, intent(inout) :: Nat      !! parent class
  type(Sparse_Con_Type), target, intent(in)    :: Acon     !! connect matrix
  type(Sparse_Val_Type), target, intent(in)    :: Aval     !! values matrix
  real,                       intent(out)   :: x(-Nat % pnt_grid % n_bnd_cells:&
                                                  Nat % pnt_grid % n_cells)
    !! unknown vector, the solution of the linear system
  real,                          intent(inout) :: b( Nat % pnt_grid % n_cells)
    !! right-hand side vector
  integer,                       intent(in)    :: miter    !! max iterations
  real,                          intent(in)    :: tol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  real, contiguous,  pointer :: r(:), p(:), q(:), d_inv(:)
  real                       :: alpha, beta, pq, rho, rho_old, res
  integer                    :: nc, iter
!==============================================================================!

  ! Take aliases
  d_inv => Aval % d_inv
  p     => Nat % p
  q     => Nat % q
  r     => Nat % r
  Grid  => Nat % pnt_grid
  nc    =  Grid % n_cells

  !----------------!
  !   r = b - Ax   !     =-->  (q used for temporary storing Ax)
  !----------------!
  call Linalg % Mat_X_Vec(nc, q, Acon, Aval, x(1:nc))  ! Ax = A * x
  call Linalg % Vec_P_Sca_X_Vec(nc, r, b, -1.0, q)     ! r  = b - Ax

  ! Check first residual
  call Linalg % Vec_D_Vec(nc, res, r, r)  ! res = r * r
  if(res < tol) return

  !---------------!
  !   z = r \ M   !    =--> (q used for z)
  !---------------!
  call Linalg % Vec_X_Vec(nc, q, r, d_inv)

  !-----------------!
  !   rho = r * z   !  =--> (q used for z)
  !-----------------!
  call Linalg % Vec_D_Vec(nc, rho, r, q)  ! rho = r * q

  !-----------!
  !   p = z   !  =--> (q used for z)
  !-----------!
  call Linalg % Vec_Copy(nc, p, q)

  do iter = 1, miter

    !------------!
    !   q = Ap   !
    !------------!
    call Linalg % Mat_X_Vec(nc, q, Acon, Aval, p)   ! q  = A * p

    !---------------------------!
    !   alfa =  rho / (p * q)   !
    !---------------------------!
    call Linalg % Vec_D_Vec(nc, pq, p, q)  ! pq = p * q
    alpha = rho / pq

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa q    !
    !---------------------!
    call Linalg % Vec_P_Sca_X_Vec(nc, x(1:nc), x(1:nc), +alpha, p)  ! x = x + alpha p
    call Linalg % Vec_P_Sca_X_Vec(nc, r,       r,       -alpha, q)  ! r = r - alpha q

    !---------------!
    !   z = r \ M   !    =--> (q used for z)
    !---------------!
    call Linalg % Vec_X_Vec(nc, q, r, d_inv)

    !-----------------!
    !   rho = r * z   !  =--> (q used for z)
    !-----------------!
    rho_old = rho
    call Linalg % Vec_D_Vec(nc, rho, r, q)  ! rho = r * q

    !--------------------------!
    !   beta = rho / rho_old   !
    !   p = z + beta * p       !  =--> (q used for p)
    !--------------------------!
    beta = rho / rho_old
    call Linalg % Vec_P_Sca_X_Vec(nc, p, q, beta, p)   ! p = q + beta p

    !--------------------!
    !   Check residual   !
    !--------------------!
    call Linalg % Vec_D_Vec(nc, res, r, r)  ! res = r * r

    if(mod(iter,32) .eq. 0) print '(a,i12,es12.3)', ' iter, res = ', iter, res
    if(res .lt. tol) goto 1

    rho_old = rho
  end do
  iter = iter - 1

1 continue
  print '(a,i12,es12.3)', ' iter, res = ', iter, res

  end subroutine
