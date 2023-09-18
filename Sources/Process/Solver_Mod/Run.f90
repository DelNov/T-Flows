!==============================================================================!
  subroutine Run(Sol, A, phi, b, norm)
!------------------------------------------------------------------------------!
!   From this procedure, the code branches either to Native or Petsc solver    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  type(Matrix_Type)  :: A
  type(Var_Type)     :: phi
  real               :: b(Sol % Nat % pnt_grid % n_cells)
  real, optional     :: norm                              ! normalization
!==============================================================================!

  ! Call linear solver to solve the equations
  if(Sol % solvers == NATIVE) then
    call Sol % Nat % Solve_Native(phi % solver,   &  ! solver
                                  phi % prec,     &  ! preconditioner
                                  A,              &
                                  phi % n,        &
                                  b,              &
                                  phi % miter,    &  ! maximum and performed...
                                  phi % niter,    &  ! ... number of iterations
                                  phi % tol,      &  ! desired tolerance
                                  phi % res,      &  ! final (achieved) resid.
                                  norm)
  else
    call Sol % Pet % Solve_Petsc(phi % solver,     &  ! solver
                                 phi % prec,       &  ! preconditioner
                                 phi % prec_opts,  &
                                 A,                &
                                 phi % n,          &
                                 b,                &
                                 phi % miter,      &  ! maximum and performed...
                                 phi % niter,      &  ! ... number of iterations
                                 phi % tol,        &  ! desired tolerance
                                 phi % res)           ! final (achieved) resid.
  end if

  end subroutine
