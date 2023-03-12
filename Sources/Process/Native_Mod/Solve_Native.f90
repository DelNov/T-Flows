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

  ! Call the desired linear solver to solve the equations
  if(solver .eq. 'bicg') then
    call Nat % Bicg(A, x, b, prec, miter, niter, tol, fin_res, norm)
  else if(solver .eq. 'cg') then
    call Nat % Cg  (A, x, b, prec, miter, niter, tol, fin_res, norm)
  else
    print *, '# ERROR: Unknown native solver: ', solver
    print *, '# This error is critical, stopping!'
    call Comm_Mod_End
    stop
  end if

  end subroutine
