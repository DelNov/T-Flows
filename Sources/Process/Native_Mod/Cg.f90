!==============================================================================!
  subroutine Cg(Nat, A, x, b, prec, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
!>  The Cg subroutine implements the Conjugate Gradient (CG) method for solving
!>  linear systems in both sequential and parallel computing environments.
!>  It solves the system Ax = b, where 'A' is a matrix, 'x' is the unknown
!>  vector, and 'b' is the right-hand side vector.  The CG method is iterative
!>  and continues until the solution converges within a specified tolerance or
!>  the maximum number of iterations (miter) is reached.  This subroutine also
!>  supports optional normalization of the system (via 'norm'), which can
!>  enhance the numerical stability and convergence behavior. Preconditioning
!>  is applied to improve convergence, with the type specified by 'prec'.
!>  The subroutine updates 'niter' with the actual number of iterations
!>  performed and 'fin_res' with the final residual value.  To avoid memory
!>  allocation and de-allocation, it relies on memory allocated in Work_Mod.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target, intent(in)    :: Nat      !! parent class
  type(Matrix_Type),  target, intent(in)    :: A        !! system matrix
  real,                       intent(out)   :: x(-Nat % pnt_grid % n_bnd_cells:&
                                                  Nat % pnt_grid % n_cells)
    !! unknown vector, the solution of the linear system
  real,                       intent(inout) :: b( Nat % pnt_grid % n_cells)
    !! right-hand side vector
  character(SL),              intent(in)    :: prec     !! preconditioner
  integer,                    intent(in)    :: miter    !! maximum iterations
  integer,                    intent(out)   :: niter    !! performed iterations
  real,                       intent(in)    :: tol      !! solver tolerance
  real,                       intent(out)   :: fin_res  !! achieved residual
  real,             optional, intent(in)    :: norm     !! normalization factor
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  integer                      :: nt, ni, nb
  real                         :: alfa, beta, rho, rho_old, bnrm2, res
  integer                      :: i, j, k, iter
  real                         :: sum_a, fn
  integer                      :: sum_n
  real,    contiguous, pointer :: p1(:), q1(:), r1(:), d(:), d_inv(:)
  real,    contiguous, pointer :: a_val(:)
  integer, contiguous, pointer :: a_col(:), a_row(:), a_dia(:)
!==============================================================================!

  call Work % Connect_Real_Cell(p1, q1, r1, d, d_inv)

  ! Take some aliases
  Grid  => Nat % pnt_grid
  a_val => A % val
  a_col => A % col
  a_row => A % row
  a_dia => A % dia

  nt = Grid % n_cells
  ni = Grid % n_cells - Grid % Comm % n_buff_cells
  nb = Grid % n_bnd_cells

  res = 0.0

  !--------------------------!
  !   Normalize the system   !
  !--------------------------!
  sum_a = 0.0
  sum_n = 0
  !$omp parallel do private(i) shared(a_val, a_dia)  &
  !$omp reduction(+ : sum_a) reduction(+ : sum_n)
  do i = 1, ni
    sum_a = sum_a + a_val(a_dia(i))
    sum_n = sum_n + 1
  end do
  !$omp end parallel do

  call Global % Sum_Real(sum_a)
  call Global % Sum_Int (sum_n)  ! this is stored somewhere, check

  sum_a = sum_a / sum_n
  fn = 1.0 / sum_a
  !$omp parallel do private(i, j) shared(a_row, a_val, b, fn)
  do i = 1, nt
    do j = a_row(i), a_row(i+1)-1
      a_val(j) = a_val(j) * fn
    end do
    b(i) = b(i) * fn
  end do
  !$omp end parallel do

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Nat % Prec_Form(ni, A, d(1:nt), d_inv(1:nt), prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Nat % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt))
  else
    bnrm2 = Nat % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Nat % Residual_Vector(ni, r1(1:nt), b(1:nt), A, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-----------!
  !   p = r   !
  !-----------!
  !$omp parallel do private(i) shared(p1, r1)
  do i = 1, ni
    p1(i) = r1(i)
  end do
  !$omp end parallel do

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, miter

    !----------------------!
    !     solve Mz = r     !
    !   (q instead of z)   !
    !----------------------!
    call Nat % Prec_Solve(ni, A, d(1:nt), d_inv(1:nt), q1(1:nt), r1(1:nt), prec)

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = 0.0
    !$omp parallel do private(i) shared(q1, r1) reduction(+ : rho)
    do i = 1, ni
      rho = rho + q1(i) * r1(i)
    end do
    !$omp end parallel do
    call Global % Sum_Real(rho)

    if(iter .eq. 1) then
      !$omp parallel do private(i) shared(p1, q1)
      do i = 1, ni
        p1(i) = q1(i)
      end do
      !$omp end parallel do
    else
      beta = rho / rho_old
      !$omp parallel do private(i) shared(p1, q1, beta)
      do i = 1, ni
        p1(i) = q1(i) + beta * p1(i)
      end do
      !$omp end parallel do
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Grid % Exchange_Cells_Real(p1(-nb:ni))
    !$omp parallel do private(i, j, k) shared(p1, q1, a_row, a_col, a_val)
    do i = 1, ni
      q1(i) = 0.0
      do j = a_row(i), a_row(i+1)-1
        k = a_col(j)
        q1(i) = q1(i) + a_val(j) * p1(k)
      end do
    end do
    !$omp end parallel do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = 0.0
    !$omp parallel do private(i) shared(p1, q1) reduction(+ : alfa)
    do i = 1, ni
      alfa = alfa + p1(i) * q1(i)
    end do
    !$omp end parallel do

    call Global % Sum_Real(alfa)
    alfa = rho / alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    !$omp parallel do private(i) shared(x, r1, p1, q1, alfa)
    do i = 1, ni
      x (i) = x (i) + alfa * p1(i)
      r1(i) = r1(i) - alfa * q1(i)
    end do
    !$omp end parallel do

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))
    else
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old = rho

  end do  ! iter

  ! Correct the number of executed iterations, because
  ! Fortran goes one number up do loop's upper limit
  iter = iter - 1

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue

  !-------------------------------------------!
  !   Refresh the solution vector's buffers   !
  !-------------------------------------------!
  call Grid % Exchange_Cells_Real(x(-nb:ni))

  !-----------------------------!
  !   De-normalize the system   !
  !-----------------------------!
  do i = 1, nt
    do j = a_row(i), a_row(i+1)-1
      a_val(j) = a_val(j) / fn
    end do
    b(i) = b(i) / fn
  end do

  fin_res = res
  niter   = iter

  call Work % Disconnect_Real_Cell(p1, q1, r1, d, d_inv)

  end subroutine
