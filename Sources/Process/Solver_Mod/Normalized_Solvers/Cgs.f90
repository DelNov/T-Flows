!==============================================================================!
  subroutine Cgs(Sol, A, x, b, prec, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a precond. CGS Method.           !
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
  use Work_Mod, only: p1         => r_cell_01,  &
                      p2         => r_cell_02,  &
                      q1         => r_cell_03,  &
                      q2         => r_cell_04,  &
                      r1         => r_cell_06,  &
                      r2         => r_cell_07,  &
                      u1         => r_cell_08,  &
                      u2         => r_cell_09,  &
                      v2         => r_cell_10,  &
                      u1_plus_q1 => r_cell_11,  &
                      fn         => r_cell_12
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined, but this         !
!   procedure is never called, so it doesn't make much sense to do it.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type), target :: Sol
  type(Matrix_Type)          :: A
  real                       :: x(-Sol % pnt_grid % n_bnd_cells :  &
                                   Sol % pnt_grid % n_cells)
  real                       :: b( Sol % pnt_grid % n_cells)
  character(SL)              :: prec     ! preconditioner
  integer                    :: miter    ! maximum and actual ...
  integer                    :: niter    ! ... number of iterations
  real                       :: tol      ! tolerance
  real                       :: fin_res  ! final residual
  real, optional             :: norm     ! normalization
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: D
  integer                    :: nt, ni, nb
  real                       :: alfa, beta, rho, rho_old, bnrm2, res
  integer                    :: i, j, k, iter
!==============================================================================!

  ! Take some aliases
  D => Sol % D
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
  call Sol % Prec_Form(ni, A, D, prec)

  !-----------------------------------!
  !    This is quite tricky point.    !
  !   What if bnrm2 is very small ?   !
  !-----------------------------------!
  if(.not. present(norm)) then
    bnrm2 = Sol % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt))
  else
    bnrm2 = Sol % Normalized_Root_Mean_Square(ni, b(1:nt), A, x(1:nt), norm)
  end if

  if(bnrm2 < tol) then
    iter=0
    goto 1
  end if

  !-----------------!
  !   r1 = b - Ax   !
  !-----------------!
  call Sol % Residual_Vector(ni, r1(1:nt), b(1:nt), A, x(1:nt))

  !--------------------------------!
  !   Calculate initial residual   !
  !--------------------------------!
  res = Sol % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))

  if(res < tol) then
    iter=0
    goto 1
  end if

  !-------------!
  !   r2 = r1   !
  !-------------!
  r2(1:ni) = r1(1:ni)

  !---------------!
  !               !
  !   Main loop   !
  !               !
  !---------------!
  do iter = 1, miter

    !-------------------!
    !   rho = (r2,z1)   !
    !-------------------!
    rho = dot_product(r1(1:ni), r2(1:ni))
    call Comm_Mod_Global_Sum_Real(rho)

    if(iter .eq. 1) then
      u1(1:ni) = r1(1:ni)
      u2(1:ni) = u1(1:ni)
    else
      beta = rho / rho_old
      u1(1:ni) = r1(1:ni) + beta *  q1(1:ni)
      u2(1:ni) = u1(1:ni) + beta * (q1(1:ni) + beta*u2(1:ni))
    end if

    !---------------------!
    !   Solve M p2 = u2   !
    !---------------------!
    call Sol % Prec_Solve(ni, A, D, p2(1:nt), u2(1:nt), prec)

    !--------------!
    !   v2 = Ap2   !
    !--------------!
    call A % pnt_grid % Exchange_Cells_Real(p2(-nb:ni))
    do i = 1, ni
      v2(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        v2(i) = v2(i) + A % val(j) * p2(k)
      end do
      alfa = alfa + r2(i) * v2(i)
    end do

    !------------------------!
    !   alfa = rho/(r2,v2)   !
    !------------------------!
    alfa = dot_product(r2(1:ni), v2(1:ni))
    call Comm_Mod_Global_Sum_Real(alfa)
    alfa=rho / alfa

    !-------------------------!
    !   q1 = u1 - alfa * v2   !
    !-------------------------!
    q1(1:ni) = u1(1:ni) - alfa * v2(1:ni)

    !-------------------------!
    !   solve Mp1 = u1 + q1   !
    !-------------------------!
    u1_plus_q1(1:ni) = u1(1:ni) + q1(1:ni)
    call Sol % Prec_Solve(ni, A, D, p1(1:nt), u1_plus_q1(1:nt), prec)

    !---------------------!
    !   x = x + alfa p1   !
    !---------------------!
    x(1:ni) = x(1:ni) + alfa * p1(1:ni)

    !---------------!
    !   q2 = A p1   !
    !---------------!
    call A % pnt_grid % Exchange_Cells_Real(p1(-nb:ni))
    do i = 1, ni
      q2(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        q2(i) = q2(i) + A % val(j) * p1(k)
      end do
    end do

    !---------------------!
    !   r = r - alfa q2   !
    !---------------------!
    r1(1:ni) = r1(1:ni) - alfa * q2(1:ni)

    !-----------------------!
    !   Check convergence   !
    !-----------------------!
    if(.not. present(norm)) then
      res = Sol % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))
    else
      res = Sol % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old=rho

  end do                ! iter

  !----------------------------------!
  !                                  !
  !   Convergence has been reached   !
  !                                  !
  !----------------------------------!
1 continue

  !-------------------------------------------!
  !   Refresh the solution vector's buffers   !
  !-------------------------------------------!
  call A % pnt_grid % Exchange_Cells_Real(x(-nb:ni))

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
