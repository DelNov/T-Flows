!==============================================================================!
  subroutine Cg(Nat, Acon, Aval, x, b, miter, niter, tol, fin_res)
!------------------------------------------------------------------------------!
!   Note: This algorithm is based on: "Templates for the Solution of Linear    !
!         Systems: Building Blocks for Iterative Methods", available for       !
!         download here: https://netlib.org/linalg/html_templates/report.html  !
!         A version which avoids "if" statement inside the Cg loop is also     !
!         implemented in the sister source, and is slighlty faster than this.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),    target, intent(inout) :: Nat      !! parent class
  type(Sparse_Con_Type), target, intent(in)    :: Acon     !! connect matrix
  type(Sparse_Val_Type), target, intent(in)    :: Aval     !! values matrix
  real, intent(out)   :: x(-Nat % pnt_grid % n_bnd_cells:&
                            Nat % pnt_grid % n_cells)
                         !! unknown vector, the solution of the linear system
  real, intent(inout) :: b( Nat % pnt_grid % n_cells)
                         !! right-hand side vector
  integer, intent(in)    :: miter    !! maximum iterations
  integer, intent(out)   :: niter    !! performed iterations
  real,    intent(in)    :: tol      !! target solver tolerance
  real,    intent(out)   :: fin_res  !! achieved residual
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  real, contiguous,  pointer :: r(:), p(:), q(:), d_inv(:)
  real                       :: fn, alpha, beta, pq, rho, rho_old, res
  integer                    :: nt, ni, iter
!==============================================================================!

  ! Take aliases
  d_inv => Aval % d_inv
  p     => Nat % p
  q     => Nat % q
  r     => Nat % r
  Grid  => Nat % pnt_grid
  nt    =  Grid % n_cells
  ni    =  Grid % n_cells - Grid % Comm % n_buff_cells

  !---------------------------------!
  !   Normalize the linear system   !
  !---------------------------------!
  call Linalg % Sys_Normalize(ni, fn, Acon, Aval, b)

  !---------------------!
  !   Preconditioning   !
  !---------------------!

  ! Scalar over diagonal (to take the mystery out: computes d_inv)
  call Linalg % Sca_O_Dia(ni, d_inv, 1.0, Acon, Aval)

  !----------------!
  !   r = b - Ax   !     =-->  (q used for temporary storing Ax)
  !----------------!
  call Linalg % Mat_X_Vec(nt, q, Acon, Aval, x(1:nt))  ! q = A * x
  call Linalg % Vec_P_Sca_X_Vec(ni, r(1:ni), b(1:ni), -1.0, q(1:ni))

  ! Check first residual
  call Linalg % Vec_D_Vec(ni, res, r(1:ni), r(1:ni))  ! res = r * r

  fin_res = res
  niter   = 0

  if(res .lt. tol) goto 2

  !-----------!
  !   p = r   !
  !-----------!
  call Linalg % Vec_Copy(ni, p(1:ni), r(1:ni))

  !--------------------------!
  !                          !
  !   Start the iterations   !
  !                          !
  !--------------------------!
  do iter = 1, miter

    !---------------!
    !   z = r \ M   !    =--> (q used for z)
    !---------------!
    call Linalg % Vec_X_Vec(ni, q(1:ni), r(1:ni), d_inv(1:ni))

    !-----------------!
    !   rho = r * z   !  =--> (q used for z)
    !-----------------!
    call Linalg % Vec_D_Vec(ni, rho, r(1:ni), q(1:ni))  ! rho = r * q

    if(iter .eq. 1) then

      !-----------!
      !   p = z   !  =--> (q used for z)
      !-----------!
      call Linalg % Vec_Copy(ni, p(1:ni), q(1:ni))  ! p = q
    else

      !--------------------------!
      !   beta = rho / rho_old   !
      !   p = z + beta * p       !  =--> (q used for p)
      !--------------------------!
      beta = rho / rho_old
      call Linalg % Vec_P_Sca_X_Vec(ni, p(1:ni), q(1:ni), beta, p(1:ni))
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Linalg % Mat_X_Vec(nt, q(1:nt), Acon, Aval, p(1:nt))

    !---------------------------!
    !   alfa =  rho / (p * q)   !
    !---------------------------!
    call Linalg % Vec_D_Vec(ni, pq, p(1:ni), q(1:ni))  ! pq = p * q
    alpha = rho / pq

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa q    !
    !---------------------!
    call Linalg % Vec_P_Sca_X_Vec(ni, x(1:ni), x(1:ni), +alpha, p(1:ni))
    call Linalg % Vec_P_Sca_X_Vec(ni, r(1:ni), r(1:ni), -alpha, q(1:ni))

    !--------------------!
    !   Check residual   !
    !--------------------!
    call Linalg % Vec_D_Vec(ni, res, r(1:ni), r(1:ni))

    if(res .lt. tol) goto 1

    rho_old = rho
  end do
  iter = iter - 1

  !-----------------------!
  !                       !
  !   End of iterations   !
  !                       !
  !-----------------------!
1 continue

  fin_res = res
  niter   = iter

  !-------------------------------------------!
  !   Refresh the solution vector's buffers   !
  !-------------------------------------------!
  call Grid % Exchange_Inside_Cells_Real(x(1:nt))

  !-------------------------------!
  !   Restore the linear system   !
  !-------------------------------!
2 continue
  call Linalg % Sys_Restore(ni, fn, Acon, Aval, b)

  end subroutine
