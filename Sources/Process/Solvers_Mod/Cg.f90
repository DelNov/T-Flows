!==============================================================================!
  subroutine Cg(mat_a, x, r1, prec, niter, tol, ini_res, fin_res, norm) 
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CG Method.            !
!------------------------------------------------------------------------------!
!   Allows preconditioning of the system by:                                   !
!     1. Diagonal preconditioning                                              !
!     2. Incomplete Cholesky preconditioning                                   !
!                                                                              !
!   The type of precondtioning is chosen by setting the variable prec to 0     !
!   (for no preconditioning), 1 (for diagonal preconditioning) or 2 (for       !
!   incomplete Cholesky preconditioning)                                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
  use Work_Mod, only: p1 => r_cell_01,  &
                      q1 => r_cell_02
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Matrix_Type) :: mat_a
  real              :: x(-mat_a % pnt_grid % n_bnd_cells :  &
                          mat_a % pnt_grid % n_cells)
  real              :: r1(mat_a % pnt_grid % n_cells)    ! [A]{x}={r1}
  character(len=80) :: prec                              ! preconditioner
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res, fin_res                  ! residual
  real, optional    :: norm                              ! normalization
!-----------------------------------[Locals]-----------------------------------!
  integer :: nt, ni, nb
  real    :: alfa, beta, rho, rho_old, bnrm2, error
  integer :: i, j, k, iter, sub
!==============================================================================!

  error = 0.

  nt = mat_a % pnt_grid % n_cells
  ni = mat_a % pnt_grid % n_cells - mat_a % pnt_grid % comm % n_buff_cells
  nb = mat_a % pnt_grid % n_bnd_cells

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(mat_a, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Normalized_Residual(ni, mat_a, x(1:nt), r1(1:ni))
  else
    bnrm2 = Normalized_Residual(ni, mat_a, x(1:nt), r1(1:ni), norm)
  end if

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual_Vector(ni, mat_a, x(1:nt), r1(1:ni))

  !-----------!
  !   p = r   !
  !-----------!
  do i = 1, ni
    p1(i) = r1(i)
  end do

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  error = Normalized_Residual(ni, mat_a, x(1:nt), r1(1:ni))

  !---------------------------------------------------------------!
  !   Residual after the correction and before the new solution   !
  !---------------------------------------------------------------!
  ini_res = error

  if(error < tol) then
    iter = 0
    goto 1
  end if

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, niter

    !----------------------!
    !     solve Mz = r     !
    !   (q instead of z)   !
    !----------------------!
    call Prec_Solve(mat_a, q1, r1, prec)

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = 0.
    do i = 1, ni
      rho = rho + r1(i)*q1(i)
    end do
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      do i = 1, ni
        p1(i) = q1(i)
      end do
    else
      beta = rho/rho_old
      do i = 1, ni
        p1(i) = q1(i) + beta*p1(i)
      end do
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Comm_Mod_Exchange_Real(mat_a % pnt_grid, p1)
    do i = 1, ni
      q1(i) = 0.
      do j = mat_a % row(i), mat_a % row(i+1)-1
        k = mat_a % col(j)
        q1(i) = q1(i) + mat_a % val(j) * p1(k)
      end do
    end do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = 0.
    do i = 1, ni
      alfa = alfa + p1(i)*q1(i)
    end do
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho/alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    do i = 1, ni
      x(i)  = x(i)  + alfa*p1(i)
      r1(i) = r1(i) - alfa*q1(i)
    end do

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      error = Normalized_Residual(ni, mat_a, x(1:nt), r1(1:ni))
    else
      error = Normalized_Residual(ni, mat_a, x(1:nt), r1(1:ni), norm)
    end if

    if(error < tol) goto 1

    rho_old = rho

  end do ! iter

1 fin_res = error
  niter = iter

  end subroutine
