!==============================================================================!
  subroutine Cg(Nat, x, miter, niter, tol, fin_res)
!------------------------------------------------------------------------------!
!   Note: This algorithm is based on: "Templates for the Solution of Linear    !
!         Systems: Building Blocks for Iterative Methods", available for       !
!         download here: https://netlib.org/linalg/html_templates/report.html  !
!         A version which avoids "if" statement inside the Cg loop is also     !
!         implemented in the sister source, and is slighlty faster than this.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target, intent(inout) :: Nat       !! parent class
  real,           intent(out) :: x(1:Nat % pnt_grid % n_cells)  !! unknown vec.
  integer,        intent(in)  :: miter    !! maximum iterations
  integer,        intent(out) :: niter    !! performed iterations
  real,           intent(in)  :: tol      !! target solver tolerance
  real,           intent(out) :: fin_res  !! achieved residual
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Sparse_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: r(:), p(:), q(:), d_inv(:), r_d_inv(:)
  real                       :: fn, alpha, beta, pq, rho, rho_old, res, bnrm2
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

  res = 0.0

  !---------------------------------!
  !   Normalize the linear system   !
  !---------------------------------!
  call Linalg % Sys_Normalize(ni, fn, A, b)

  !---------------------!
  !   Preconditioning   !
  !---------------------!

  ! Scalar over diagonal (to take the mystery out: computes d_inv)
  call Linalg % Sca_O_Dia(ni, d_inv, 1.0, A)

  bnrm2 = Linalg % Normalized_Root_Mean_Square(ni, b(1:nt), A)

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !     =-->  (q used for temporary storing Ax)
  !----------------!
  call Linalg % Mat_X_Vec(nt, q(1:ni), A, x(1:nt))  ! q = A * x
  call Linalg % Vec_P_Sca_X_Vec(ni, r(1:ni), b(1:ni), -1.0, q(1:ni))

  ! Check first residual
  res = Linalg % Normalized_Root_Mean_Square(ni, r(1:nt), A)

  if(res < tol) then
    iter = 0
    goto 1
  end if

  fin_res = res
  niter   = 0

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

    !--------------------!
    !   Check residual   !
    !--------------------!
    res = Linalg % Normalized_Root_Mean_Square(ni, r(1:nt), A)

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
  call Linalg % Sys_Restore(ni, fn, A, b)

  call Work % Disconnect_Real_Cell(p, q, r, r_d_inv)

  end subroutine
