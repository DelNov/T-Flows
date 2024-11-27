!==============================================================================!
  subroutine Prec_Solve(Nat, ni, A, d, d_inv, x, b, prec)
!------------------------------------------------------------------------------!
!>  The Prec_Solve subroutine is integral for solving preconditioned systems
!>  in iterative solvers. It applies the preconditioning matrix 'd' to the
!>  system [d]{x}={b}, aiding in improving solver efficiency and convergence.
!>  The subroutine caters to different preconditioning strategies: Diagonal
!>  (Jacobi) preconditioning, Incomplete Cholesky preconditioning, or no
!>  preconditioning. The choice of the preconditioning method is determined by
!>  the 'prec' parameter. For Diagonal preconditioning, it solves the system
!>  element-wise. Incomplete Cholesky involves a more complex calculation with
!>  forward and backward substitution.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),         intent(in)  :: Nat       !! parent class
  integer,                    intent(in)  :: ni        !! number of unknowns
  type(Matrix_Type),  target, intent(in)  :: A         !! system matrix
  real,                       intent(in)  :: d(:)      !! preconditioned matrix
  real,                       intent(in)  :: d_inv(:)  !! preconditioned inverse
  real,                       intent(out) :: x(:)      !! unknown vector
  real,                       intent(in)  :: b(:)      !! right-hand side vector
  character(SL),              intent(in)  :: prec      !! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  integer                       :: i, j, k
  real                          :: sum
  integer, contiguous,  pointer :: a_col(:), a_row(:), a_dia(:)
  real,    contiguous,  pointer :: a_val(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
!==============================================================================!

  call Profiler % Start('Native_Prec_Solve (all solvers)')

  ! Take some aliases
  a_col => A % col
  a_row => A % row
  a_dia => A % dia
  a_val => A % val

  !---------------------------------!
  !   1) diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'jacobi') then
    !$omp parallel do private(i) shared (x, b, d_inv)
    do i = 1, ni
      x(i) = b(i) * d_inv(i)
    end do
    !$omp end parallel do

  !--------------------------------------------!
  !   2) incomplete cholesky preconditioning   !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then

    ! Forward substitutionn
    do i = 1, ni
      sum = b(i)
      do j = a_row(i), a_dia(i) - 1              ! only the lower triangular
        k = a_col(j)
        sum = sum - a_val(j) * x(k)
      end do
      x(i) = sum * d(i)
    end do

    !$omp parallel do private(i) shared (x, d_inv)
    do i = 1, ni
      x(i) = x(i) * d_inv(i)
    end do
    !$omp end parallel do

    ! Backward substitution
    do i = ni, 1, -1
      sum = x(i)
      do j = a_dia(i) + 1, a_row(i+1) - 1        ! upper triangular
        k = a_col(j)
        if(k <= ni) sum = sum - a_val(j) * x(k)  ! avoid buffer entries
      end do
      x(i) = sum * d(i)
    end do

  !---------------------------!
  !   .) no preconditioning   !
  !---------------------------!
  else
    !$omp parallel do private(i) shared (x, b)
    do i = 1, ni
      x(i) = b(i)
    end do
    !$omp end parallel do
  end if

  call Profiler % Stop('Native_Prec_Solve (all solvers)')

  end subroutine
