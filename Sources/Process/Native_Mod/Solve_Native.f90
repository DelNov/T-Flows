!==============================================================================!
  subroutine Solve_Native(Nat,           &
                          solver, prec,  &  ! strings
                          A, x, b,       &  ! the whole system
                          miter, niter,  &  ! iterations
                          tol, fin_res,  &  ! tolerance, residual
                          norm)             ! normalization
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Nat
  character(*)       :: solver   ! solver
  character(*)       :: prec     ! preconditioner
  type(Matrix_Type)  :: A
  real               :: x(-Nat % pnt_grid % n_bnd_cells :  &
                           Nat % pnt_grid % n_cells)
  real               :: b( Nat % pnt_grid % n_cells)
  integer            :: miter    ! maximum and actual ...
  integer            :: niter    ! ... number of iterations
  real               :: tol      ! tolerance
  real               :: fin_res  ! final residual
  real, optional     :: norm     ! normalization
!==============================================================================!

  PRINT *, '@SOLVE_NATIVE: ', solver, ' ', prec

  ! Call the desired linear solver to solve the equations
  if(solver .eq. 'BICG') then
    call Nat % Bicg(A, x, b, prec, miter, niter, tol, fin_res, norm)
  else if(solver .eq. 'CG') then
    call Nat % Cg  (A, x, b, prec, miter, niter, tol, fin_res, norm)
  else if(solver .eq. 'CGS') then
    call Nat % Cgs (A, x, b, prec, miter, niter, tol, fin_res, norm)
  end if

  end subroutine
