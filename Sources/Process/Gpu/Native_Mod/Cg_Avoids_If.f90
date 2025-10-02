!==============================================================================!
  subroutine Cg(Nat, x, miter, niter, tol, fin_res)
!------------------------------------------------------------------------------!
!   Note: This is an alternative algorithm and I am honestly not sure where    !
!         I have found it any more, but it avoids one "if" block during the    !
!         iterations.  It gives exactly the same result as the original Cg.    !
!         although it is slighly faster.  This verision was implemented in     !
!         SFS, and I wanted to make sure it is the same as the original Cg.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target, intent(inout) :: Nat       !! parent class
  real, intent(out)    :: x(1:Nat % pnt_grid % n_cells)  !! unknown vector
  integer, intent(in)  :: miter    !! maximum iterations
  integer, intent(out) :: niter    !! performed iterations
  real,    intent(in)  :: tol      !! target solver tolerance
  real,    intent(out) :: fin_res  !! achieved residual
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Sparse_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: r(:), p(:), q(:), d_inv(:), r_d_inv(:)
  real                       :: fn, alpha, beta, pq, rho, rho_old, res
  integer                    :: nt, ni, iter
!==============================================================================!

  call Work % Connect_Real_Cell(p, q, r, r_d_inv)

  ! Take aliases
  Grid  => Nat % pnt_grid
  A     => Nat % A
  d_inv => A % d_inv
  b     => Nat % b
  nt    =  Grid % n_cells
  ni    =  Grid % n_cells - Grid % Comm % n_buff_cells

  !---------------------------------!
  !   Normalize the linear system   !
  !---------------------------------!
  call Linalg % Sys_Normalize(ni, fn, A, b)

  !---------------------!
  !   Preconditioning   !
  !---------------------!

  ! Scalar over diagonal (to take the mystery out: computes d_inv)
  call Linalg % Sca_O_Dia(ni, d_inv, 1.0, A)

  !----------------!
  !   r = b - Ax   !     =-->  (q used for temporary storing Ax)
  !----------------!
  call Linalg % Mat_X_Vec(nt, q(1:ni), A, x(1:nt))  ! q = A * x
  call Linalg % Vec_P_Sca_X_Vec(ni, r(1:ni), b(1:ni), -1.0, q(1:ni))

  ! Check first residual
  call Linalg % Vec_X_Vec(ni, r_d_inv(1:ni), r(1:ni), d_inv(1:ni))
  call Linalg % Vec_D_Vec(ni, res, r_d_inv(1:ni), r_d_inv(1:ni))
  res = sqrt(res)

  fin_res = res
  niter   = 0

  if(res .lt. tol) goto 2

  !---------------!
  !   z = r \ M   !    =--> (q used for z)
  !---------------!
  call Linalg % Vec_X_Vec(ni, q(1:ni), r(1:ni), d_inv(1:ni))

  !-----------------!
  !   rho = r * z   !  =--> (q used for z)
  !-----------------!
  call Linalg % Vec_D_Vec(ni, rho, r(1:ni), q(1:ni))  ! rho = r * q

  !-----------!
  !   p = z   !  =--> (q used for z)
  !-----------!
  call Linalg % Vec_Copy(ni, p(1:ni), q(1:ni))

  !--------------------------!
  !                          !
  !   Start the iterations   !
  !                          !
  !--------------------------!
  do iter = 1, miter

    !------------!
    !   q = Ap   !
    !------------!
    call Linalg % Mat_X_Vec(nt, q(1:ni), A, p(1:nt))

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

    !---------------!
    !   z = r \ M   !    =--> (q used for z)
    !---------------!
    call Linalg % Vec_X_Vec(ni, q(1:ni), r(1:ni), d_inv(1:ni))

    !-----------------!
    !   rho = r * z   !  =--> (q used for z)
    !-----------------!
    rho_old = rho
    call Linalg % Vec_D_Vec(ni, rho, r(1:ni), q(1:ni))

    !--------------------------!
    !   beta = rho / rho_old   !
    !   p = z + beta * p       !  =--> (q used for p)
    !--------------------------!
    beta = rho / rho_old
    call Linalg % Vec_P_Sca_X_Vec(ni, p(1:ni), q(1:ni), beta, p(1:ni))

    !--------------------!
    !   Check residual   !
    !--------------------!
    call Linalg % Vec_X_Vec(ni, r_d_inv(1:ni), r(1:ni), d_inv(1:ni))
    call Linalg % Vec_D_Vec(ni, res, r_d_inv(1:ni), r_d_inv(1:ni))
    res = sqrt(res)

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
  call Linalg % Sys_Restore(ni, fn, A, b)

  call Work % Disconnect_Real_Cell(p, q, r, r_d_inv)

  end subroutine
