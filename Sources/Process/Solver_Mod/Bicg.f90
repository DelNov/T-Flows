!==============================================================================!
  subroutine Bicg(sol, x, b, prec, niter, tol, ini_res, fin_res, norm) 
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. BiCG Method.          !
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
                      p2 => r_cell_02,  &
                      q1 => r_cell_03,  &
                      q2 => r_cell_04,  &
                      r1 => r_cell_05,  &
                      r2 => r_cell_06
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type), target :: sol
  real              :: x(-sol % pnt_grid % n_bnd_cells :  &
                          sol % pnt_grid % n_cells)
  real              :: b( sol % pnt_grid % n_cells)      ! [A]{x}={b}
  character(len=80) :: prec                              ! preconditioner
  integer           :: niter                             ! number of iterations
  real              :: tol                               ! tolerance
  real              :: ini_res, fin_res                  ! residual
  real, optional    :: norm                              ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: a
  type(Matrix_Type), pointer :: d
  integer                    :: nt, ni, nb
  real                       :: alfa, beta, rho, rho_old, bnrm2, res
  integer                    :: i, j, k, iter, sub
!==============================================================================!

  ! Take some aliases
  a => sol % a
  d => sol % d
  nt = a % pnt_grid % n_cells
  ni = a % pnt_grid % n_cells - a % pnt_grid % comm % n_buff_cells
  nb = a % pnt_grid % n_bnd_cells

  res = 0.0

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(ni, a, d, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Normalized_Root_Mean_Square(ni, b(1:nt), a, x(1:nt))
  else
    bnrm2 = Normalized_Root_Mean_Square(ni, b(1:nt), a, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual_Vector(ni, r1(1:nt), b(1:nt), a, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Normalized_Root_Mean_Square(ni, r1(1:nt), a, x(1:nt))
  ini_res = res

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-----------------------!
  !   Choose initial r~   !
  !-----------------------!
  do i = 1, ni
    r2(i) = r1(i)
  end do

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, niter

    !------------------------!
    !    solve M   z  = r    !
    !    solve M^T z~ = r~   !  don't have M^T!!!
    !    (q instead of z)    !
    !------------------------!
    call Prec_Solve(ni, a, d, q1(1:nt), r1(1:nt), prec)
    call Prec_Solve(ni, a, d, q2(1:nt), r2(1:nt), prec)

    !------------------!
    !   rho = (z,r~)   !
    !------------------!
    rho = dot_product(q1(1:ni), r2(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      p1(1:ni) = q1(1:ni)
      p2(1:ni) = q2(1:ni)
    else
      beta = rho / rho_old
      p1(1:ni) = q1(1:ni) + beta * p1(1:ni)
      p2(1:ni) = q2(1:ni) + beta * p2(1:ni)
    end if

    !---------------!
    !   q = A   p   !
    !   q~= A^T p~  !  don't have A^T
    !---------------!
    call Comm_Mod_Exchange_Real(a % pnt_grid, p1)
    call Comm_Mod_Exchange_Real(a % pnt_grid, p2)
    do i = 1, ni
      q1(i) = 0.0
      q2(i) = 0.0
      do j = a % row(i), a % row(i+1)-1
        k = a % col(j)
        q1(i) = q1(i) + a % val(j) * p1(k)
        q2(i) = q2(i) + a % val(j) * p2(k)
      end do
    end do

    !----------------------!
    !   alfa = rho/(p,q)   !
    !----------------------!
    alfa = 0.0
    alfa = alfa + dot_product(p2(1:ni), q1(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho / alfa

    !--------------------!
    !   x = x + alfa p   !
    !   r = r - alfa q   !
    !--------------------!
    x (1:ni) = x (1:ni) + alfa*p1(1:ni)
    r1(1:ni) = r1(1:ni) - alfa*q1(1:ni)
    r2(1:ni) = r2(1:ni) - alfa*q2(1:ni)

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      res = Normalized_Root_Mean_Square(ni, r1(1:nt), a, x(1:nt))
    else
      res = Normalized_Root_Mean_Square(ni, r1(1:nt), a, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old=rho

  end do     ! iter

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue
  fin_res = res
  niter = iter

  end subroutine
