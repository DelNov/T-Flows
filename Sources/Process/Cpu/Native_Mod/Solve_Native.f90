!==============================================================================!
  subroutine Solve_Native(Nat,           &
                          solver, prec,  &  ! strings
                          A, x, b,       &  ! the whole system
                          miter, niter,  &  ! iterations
                          tol, fin_res)     ! tolerance, residual
!------------------------------------------------------------------------------!
!>  Solve_Native is a decision-making subroutine in the Native_Mod module that
!>  determines which native linear solver to use based on user preferences.
!>  It acts as a branch point, delegating the task of solving the linear
!>  systems of equations to the appropriate solver.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), intent(in)    :: Nat      !! parent class
  character(*),       intent(in)    :: solver   !! solver to use
  character(*),       intent(in)    :: prec     !! preconditioner to use
  type(Matrix_Type),  intent(in)    :: A        !! system matrix
  real,               intent(out)   :: x(-Nat % pnt_grid % n_bnd_cells :  &
                                          Nat % pnt_grid % n_cells)
    !! unknown solution vecctor
  real,               intent(inout) :: b( Nat % pnt_grid % n_cells)
    !! right hand side vector
  integer,            intent(in)    :: miter    !! maximum iterations
  integer,            intent(out)   :: niter    !! performed
  real,               intent(in)    :: tol      !! target tolerance
  real,               intent(out)   :: fin_res  !! final (achieved) residual
!==============================================================================!

  ! Call the desired linear solver to solve the equations
  if(solver .eq. 'bicg') then
    call Nat % Bicg(A, x, b, prec, miter, niter, tol, fin_res)
  else if(solver .eq. 'cg') then
    call Nat % Cg  (A, x, b, prec, miter, niter, tol, fin_res)
  else
    call Message % Error(64, 'Unknown native solver: '//trim(solver)  //  &
                             '. This error is critical, exiting! '    //  &
                             'Check the file: Documents/all_control'  //  &
                             '_keywords to see which solvers are '    //  &
                             'available.',                                &
                             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
