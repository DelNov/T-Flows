!==============================================================================!
  subroutine Run(Sol, solver, prec, A, x, b, miter, niter, tol, fin_res, norm)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  character(*)       :: solver   ! solver
  character(*)       :: prec     ! preconditioner
  type(Matrix_Type)  :: A
  real               :: x(-Sol % Nat % pnt_grid % n_bnd_cells :  &
                           Sol % Nat % pnt_grid % n_cells)
  real               :: b( Sol % Nat % pnt_grid % n_cells)
  integer            :: miter    ! maximum and actual ...
  integer            :: niter    ! ... number of iterations
  real               :: tol      ! tolerance
  real               :: fin_res  ! final residual
  real, optional     :: norm     ! normalization
!==============================================================================!

  ! Call linear solver to solve the equations
  if(Sol % solvers == NATIVE) then
    if(solver .eq. 'CG') then
      call Sol % Nat % Cg_(A,        &
                          x,        &
                          b,        &
                          prec,     &
                          miter,    &
                          niter,    &
                          tol,      &
                          fin_res,  &
                          norm)
    else if(solver .eq. 'BICG') then
      call Sol % Nat % Bicg(A,        &
                            x,        &
                            b,        &
                            prec,     &
                            miter,    &
                            niter,    &
                            tol,      &
                            fin_res,  &
                            norm)
    end if
  else
    call Sol % Pet % Solve(Sol % Nat,  &
                           A,          &
                           x,          &
                           b,          &
                           miter,      &
                           niter,      &
                           tol)
  end if

  end subroutine
