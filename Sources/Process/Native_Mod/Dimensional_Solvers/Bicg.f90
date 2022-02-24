!==============================================================================!
  subroutine Bicg(Nat, A, x, b, prec, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
!   Solves the linear systems of equations by a preconditioned BiCG method     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: p1 => r_cell_11,  &
                      p2 => r_cell_12,  &
                      q1 => r_cell_13,  &
                      q2 => r_cell_14,  &
                      r1 => r_cell_15,  &
                      r2 => r_cell_16
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                  (allocates Work_Mod)                             !
!     |                                                                        !
!     +----> Compute_Energy   (uses r_cell_01..03)                             !
!     |        |                                                               !
!     +----> Compute_Momentum (does not use Work_Mod)                          !
!     |        |                                                               !
!     +----> Compute_Scalar   (uses r_cell_04)                                 !
!              |                                                               !
!              +----> Bicg    (safe to use r_cell_11..16)                      !
!                                                                              !
!   Main_Pro                                    (allocates Work_Mod)           !
!     |                                                                        !
!     +----> Turb_Mod_Main                      (does not use Work_Mod)        !
!              |                                                               !
!              +---> Turb_Mod_Compute_Variable  (does not use Work_Mod)        !
!              |       |                                                       !
!              +---> Turb_Mod_Compute_Stress    (uses r_cell_01..09)           !
!                      |                                                       !
!                      +----> Bicg              (safe to use r_cell_11..16)    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target :: Nat
  type(Matrix_Type)          :: A
  real                       :: x(-Nat % pnt_grid % n_bnd_cells :  &
                                   Nat % pnt_grid % n_cells)
  real                       :: b( Nat % pnt_grid % n_cells)
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
  D => Nat % D
  nt = A % pnt_grid % n_cells
  ni = A % pnt_grid % n_cells - A % pnt_grid % comm % n_buff_cells
  nb = A % pnt_grid % n_bnd_cells

  res = 0.0

  !---------------------!
  !   Preconditioning   !
  !---------------------!
  call Nat % Prec_Form(ni, A, D, prec)

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
  do iter = 1, miter

    !------------------------!
    !    solve M   z  = r    !
    !    solve M^T z~ = r~   !  don't have M^T!!!
    !    (q instead of z)    !
    !------------------------!
    call Nat % Prec_Solve(ni, A, D, q1(1:nt), r1(1:nt), prec)
    call Nat % Prec_Solve(ni, A, D, q2(1:nt), r2(1:nt), prec)

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
    call A % pnt_grid % Exchange_Cells_Real(p1(-nb:ni))
    call A % pnt_grid % Exchange_Cells_Real(p2(-nb:ni))
    do i = 1, ni
      q1(i) = 0.0
      q2(i) = 0.0
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        q1(i) = q1(i) + A % val(j) * p1(k)
        q2(i) = q2(i) + A % val(j) * p2(k)
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
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt))
    else
      res = Nat % Normalized_Root_Mean_Square(ni, r1(1:nt), A, x(1:nt), norm)
    end if

    if(res < tol) goto 1

    rho_old = rho

  end do ! iter

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
  call A % pnt_grid % Exchange_Cells_Real(x(-nb:ni))

  fin_res = res
  niter   = iter

  end subroutine
