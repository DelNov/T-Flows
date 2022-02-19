!==============================================================================!
  subroutine Run(Sol, solver, prec, A, x, b, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
!   From this procedure, the code branches either to Native or Petsc solver    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  character(*)       :: solver   ! solver and preconditioner as ...
  character(*)       :: prec     ! ... specifed in the control file
  type(Matrix_Type)  :: A
  real               :: x(-Sol % Nat % pnt_grid % n_bnd_cells :  &
                           Sol % Nat % pnt_grid % n_cells)
  real               :: b( Sol % Nat % pnt_grid % n_cells)
  integer            :: miter    ! maximum and actual ...
  integer            :: niter    ! ... number of iterations
  real               :: tol      ! desired tolerance
  real               :: fin_res  ! final (achieved) residual
  real, optional     :: norm     ! normalization
!==============================================================================!

  ! Call linear solver to solve the equations
  if(Sol % solvers == NATIVE) then
    call Sol % Nat % Solve_Native(solver,   &
                                  prec,     &
                                  A,        &
                                  x,        &
                                  b,        &
                                  miter,    &
                                  niter,    &
                                  tol,      &
                                  fin_res,  &
                                  norm)
  else
    call Sol % Pet % Solve_Petsc(Sol % Nat,  &  ! 1st argument native solver
                                 solver,     &
                                 prec,       &
                                 A,          &
                                 x,          &
                                 b,          &
                                 miter,      &
                                 niter,      &
                                 tol)
  end if

  end subroutine
