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
  class(Native_Type), intent(in)    :: Nat
  character(*),       intent(in)    :: solver   ! solver
  character(*),       intent(in)    :: prec     ! preconditioner
  type(Matrix_Type),  intent(in)    :: A
  real,               intent(out)   :: x(-Nat % pnt_grid % n_bnd_cells :  &
                                          Nat % pnt_grid % n_cells)
  real,               intent(inout) :: b( Nat % pnt_grid % n_cells)
  integer,            intent(in)    :: miter    ! maximum and actual ...
  integer,            intent(out)   :: niter    ! ... number of iterations
  real,               intent(in)    :: tol      ! tolerance
  real,               intent(out)   :: fin_res  ! final residual
  real,     optional, intent(in)    :: norm     ! normalization
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
