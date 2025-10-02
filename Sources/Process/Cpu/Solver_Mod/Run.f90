!==============================================================================!
  subroutine Run(Sol, A, phi, b, norm)
!------------------------------------------------------------------------------!
!>  The Run subroutine within the Solver_Mod module is the best example on how
!>  encapsulation of decision-making process for branching between native
!>  solvers and PETSc solvers works in Solver_Mod. This approach effectively
!>  abstracts away the underlying solver details from the rest of the code,
!>  providing a unified interface regardless of the chosen solver type.  In
!>  other words: this subroutine serves as a gateway to either native or PETSc
!>  solvers, based on the solver type specified in Sol % solvers.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol                                !! parent class
  type(Matrix_Type)  :: A                                  !! T-Flows matrix
  type(Var_Type)     :: phi                                !! unknown variable
  real               :: b(Sol % Nat % pnt_grid % n_cells)  !! source vector
  real, optional     :: norm                               !! normalization
                                                           !! factor
!==============================================================================!

  ! Call linear solver to solve the equations
  if(Sol % solvers == NATIVE) then
    call Sol % Nat % Solve_Native(phi % solver,  &  !! solver's name
                                  phi % prec,    &  !! preconditioner's name
                                  A,             &  !! T-Flows matrix object
                                  phi % n,       &  !! unknown vector
                                  b,             &  !! source vector
                                  phi % miter,   &  !! maximum iterations
                                  phi % niter,   &  !! performed iterations
                                  phi % tol,     &  !! target tolerance
                                  phi % res,     &  !! final/achieved residual
                                  norm)             !! normalization factor
  else
    call phi % Pet % Solve_Petsc(phi % solver,   &  !! solver's name
                                 phi % prec,     &  !! preconditioner's name
                                 phi % o_prec,   &  !! preconditioner's options
                                 A,              &  !! T-Flows matrix object
                                 phi % n,        &  !! unknown vector
                                 b,              &  !! source vector
                                 phi % miter,    &  !! maximum iterations
                                 phi % niter,    &  !! performed iterations
                                 phi % tol,      &  !! target tolerance
                                 phi % res,      &  !! final/achieved residual
                                 phi % blend_matrix)  !! flag telling if matrix
                                                      !! is upwind blended
  end if

  end subroutine
