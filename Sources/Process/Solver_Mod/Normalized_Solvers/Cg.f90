!==============================================================================!
  subroutine Solver_Mod_Cg(sol, x, b, prec, miter, niter, tol, fin_res, norm)
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
  use Work_Mod, only: p1 => r_cell_01,  &
                      q1 => r_cell_02,  &
                      r1 => r_cell_03,  &
                      fn => r_cell_04
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Solver_Type), target :: sol
  real              :: x(-sol % pnt_grid % n_bnd_cells :  &
                          sol % pnt_grid % n_cells)
  real              :: b( sol % pnt_grid % n_cells)  ! [A]{x}={b}
  character(SL)     :: prec                          ! preconditioner
  integer           :: miter                         ! max and actual ...
  integer           :: niter                         ! ... num. iterations
  real              :: tol                           ! tolerance
  real              :: fin_res                       ! residual
  real, optional    :: norm                          ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: A
  type(Matrix_Type), pointer :: d
  integer                    :: nt, ni, nb
  real                       :: alfa, beta, rho, rho_old, bnrm2, res
  integer                    :: i, j, k, iter
!==============================================================================!

  ! Take some aliases
  A => sol % A
  d => sol % d
  nt = A % pnt_grid % n_cells
  ni = A % pnt_grid % n_cells - A % pnt_grid % comm % n_buff_cells
  nb = A % pnt_grid % n_bnd_cells

  res = 0.0

  !--------------------------!
  !   Normalize the system   !
  !--------------------------!
  do i = 1, nt
    fn(i) = 1.0 / A % val(A % dia(i))
    do j = A % row(i), A % row(i+1)-1
      A % val(j) = A % val(j) * fn(i)
    end do
    b(i) = b(i) * fn(i)
  end do

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Prec_Form(ni, A, d, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt))
  else
    bnrm2 = Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then
    iter = 0
    goto 1
  end if

  !----------------!
  !   r = b - Ax   !
  !----------------!
  call Residual_Vector(ni, r1(1:nt), b(1:nt), A, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))

  if(res < tol) then
    iter = 0
    goto 1
  end if

  !-----------!
  !   p = r   !
  !-----------!
  p1(1:ni) = r1(1:ni)

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
    call Prec_Solve(ni, A, d, q1(1:nt), r1(1:nt), prec)

    !-----------------!
    !   rho = (r,z)   !
    !-----------------!
    rho = dot_product(r1(1:ni), q1(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      p1(1:ni) = q1(1:ni)
    else
      beta = rho / rho_old
      p1(1:ni) = q1(1:ni) + beta * p1(1:ni)
    end if

    !------------!
    !   q = Ap   !
    !------------!
    call Grid_Mod_Exchange_Cells_Real(A % pnt_grid, p1(-nb:ni))
    do i = 1, ni
      q1(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        q1(i) = q1(i) + A % val(j) * p1(k)
      end do
    end do

    !------------------------!
    !   alfa = (r,z)/(p,q)   !
    !------------------------!
    alfa = dot_product(p1(1:ni), q1(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa = rho/alfa

    !---------------------!
    !   x = x + alfa p    !
    !   r = r - alfa Ap   !
    !---------------------!
    x (1:ni) = x (1:ni) + alfa * p1(1:ni)
    r1(1:ni) = r1(1:ni) - alfa * q1(1:ni)

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      res = Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))
    else
      res = Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old = rho

  end do ! iter

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue

  !-----------------------------!
  !   De-normalize the system   !
  !-----------------------------!
  do i = 1, nt
    do j = A % row(i), A % row(i+1)-1
      A % val(j) = A % val(j) / fn(i)
    end do
    b(i) = b(i) / fn(i)
  end do

  fin_res = res
  niter   = iter

  end subroutine
