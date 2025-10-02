!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         blend_matrix)
!------------------------------------------------------------------------------!
!>  The subroutine Solve_Petsc in the Petsc_Mod module is an integral
!>  component for solving linear systems in T-Flows using PETSc solvers.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Matrix and vector preparation:                                           !
!     - Conditionally copies values from T-Flows' matrix A to PETSc's matrix   !
!       Pet % A. This step is necessary for updating the PETSc matrix with     !
!       the latest values from the simulation and is performed if the matrix   !
!       wasn't previously copied or if the blend_matrix flag is true.          !
!     - Fills PETSc vectors Pet % x (solution vector) and Pet % b (right-hand  !
!       side) with values from T-Flows' vectors x and b.                       !
!   * Solver and preconditioner setup:                                         !
!     - Sets up the PETSc Krylov Subspace (KSP) solver and preconditioner      !
!       based on user-defined options (solver, prec).                          !
!     - Handles the specification of preconditioner options (prec_opts),       !
!       allowing for granular control over the preconditioning process.        !
!     - Sets the initial guess for the solver to non-zero, allowing the solver !
!       to start from the existing solution estimate, which can enhance        !
!       convergence in iterative simulations.                                  !
!   * Solver execution:                                                        !
!     - Sets tolerances for the solver based on the provided tol value and     !
!       the maximum number of iterations miter.                                !
!     - Calls the PETSc solver to solve the linear system. This process        !
!       involves iterative methods and can be computationally intensive,       !
!       depending on the complexity of the system.                             !
!     - Retrieves the number of iterations performed (niter) and the final     !
!       residual norm (fin_res), providing insights into the solver's          !
!       performance and convergence.                                           !
!   * Solution retrieval:                                                      !
!     - Upon successful convergence, copies the solution from the PETSc vector !
!       back to T-Flows' vector x. This step ensures that the computed         !
!       solution is available in T-Flows for further processing or analysis.   !
!   * Efficiency considerations:                                               !
!     - Profiling is used around key sections (matrix copying, preconditioner  !
!       setup) to monitor performance and identify potential bottlenecks.      !
!     - The subroutine is optimized to minimize unnecessary operations. For    !
!       instance, matrix copying and preconditioner setup are conditionally    !
!       executed to avoid redundant computations.                              !
!   * Robustness and flexibility:                                              !
!     - The subroutine is designed to be robust against solver convergence     !
!       issues, checking the convergence reason and handling cases where the   !
!       solver does not converge.                                              !
!     - The flexible design allows for various solver and preconditioner       !
!       configurations, making it adaptable to different simulation scenarios  !
!       and requirements.                                                      !
!   * Interaction with PETSc:                                                  !
!     - Extensive use of C_Petsc_* interface functions for direct interaction  !
!       with PETSc's C API, reflecting a consistent approach to integrate      !
!       PETSc functionalities within T-Flows.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)          :: Pet      !! parent object of the Petsc_Type
  character(*),  intent(in)  :: solver   !! name of the solver to use
  character(*),  intent(in)  :: prec     !! name of the preconditioner to use
  character(SL), intent(in)  :: prec_opts(MAX_STRING_ITEMS)
    !! list of options passed to preconditioner from the control file
  type(Matrix_Type)          :: A        !! matrix in T-Flows format
  real                       :: x(-Pet % pnt_grid % n_bnd_cells :  &
                                   Pet % pnt_grid % n_cells)  !! unknown vector
  real                       :: b( Pet % pnt_grid % n_cells)  !! source vector
  integer,       intent(in)  :: miter    !! maximum number of iterations
  integer,       intent(out) :: niter    !! performed number of iterations
  real,          intent(in)  :: tol      !! target solver tolerance
  real,          intent(out) :: fin_res  !! final residual after linear solver
  logical,       intent(in)  :: blend_matrix
    !! flag indicating if the system matrix is blended with upwind terms
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, k, j_loc, a_glo_j(256)
  character(SL) :: solvers  ! fortran string to store solver
  character(SL) :: precs    ! fortran string to store preconditioner
  integer       :: l        ! length of the two strings above
!==============================================================================!

  !-----------------------------------------------------------!
  !   Fill up PETSc matrix with values from original matrix   !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !   (This can be very slow, comparable to call to solver)   !
  !-----------------------------------------------------------!
  if(.not. Pet % matrix_coppied .or. blend_matrix) then
    call Profiler % Start('Solve_Petsc (matrix copy)')

    do i = 1, Pet % m_lower
      do j = A % row(i), A % row(i+1)-1
        k = A % col(j)
        j_loc = j - A % row(i) + 1   ! local j counter for this row
        a_glo_j(j_loc) = A % glo(k)  ! fill up the column global number
      end do
      call C_Petsc_Mat_Set_Values(Pet % A,     &        ! matrix
                                  1,           &        ! one row at a time
                                  A % glo(i),  &        ! global row number
                                  j_loc,       &        ! number of columns
                                  a_glo_j,     &        ! global column numbers
                                  A % val(A % row(i)))  ! matrix entries
    end do

    ! The following two calls are needed after calls to MatSetValue
    ! But be carefull when you call: it triggers the re-formation of
    ! the preconditioning matrices which is a rather slow process
    call C_Petsc_Mat_Assemble(Pet % A)

    Pet % matrix_coppied = .true.
    call Profiler % Stop('Solve_Petsc (matrix copy)')
  end if

  !---------------------!
  !   Fill up vectors   ! (A % glo starts from zero)
  !---------------------!
  do i = 1, Pet % m_lower
    call C_Petsc_Vec_Set_Value(Pet % x, A % glo(i), x(i))
    call C_Petsc_Vec_Set_Value(Pet % b, A % glo(i), b(i))
  end do

  ! The following two calls are needed after the calls to VecSetValue
  call C_Petsc_Vec_Assemble(Pet % x)
  call C_Petsc_Vec_Assemble(Pet % b)

  !-----------------------------------!
  !   Set solver and preconditioner   !
  !- - - - - - - - - - - - - - - - - -!
  !   (This can be very slow, much    !
  !    slower than call to solver.)   !
  !-----------------------------------!
  if(.not. Pet % precond_formed .or. blend_matrix) then
    call Profiler % Start('Solve_Petsc (precondition)')
    call C_Petsc_Ksp_Set_Operators(Pet % ksp,  &  ! solver
                                   Pet % A)       ! matrix

    ! Set (choose) the Krylov subspace method to be used
    solvers=solver;  l=len_trim(solvers);  solvers(l+1:l+1)=c_null_char
    call C_Petsc_Ksp_Set_Type(Pet % ksp,  &  ! solver
                              solvers)

    ! Set (choose) the preconditioner to be used
    precs=prec;  l=len_trim(precs);  precs(l+1:l+1)=c_null_char
    call C_Petsc_Ksp_Set_Preconditioner(Pet % ksp,  &  ! solver
                                        Pet % pc,   &  ! preconditioner
                                        precs)
    Pet % precond_formed = .true.
    call Profiler % Stop('Solve_Petsc (precondition)')
  end if

  !--------------------------------------------!
  !   Do not start from zero as inital guess   !
  !--------------------------------------------!
  call C_Petsc_Ksp_Set_Initial_Guess_Nonzero(Pet % ksp)

  !------------------------------------!
  !   Process preconditioner options   !
  !------------------------------------!
  if(prec_opts(1) .ne. '') then
    i = 1
    do while(i < MAX_STRING_ITEMS .and. prec_opts(i)(1:1) .ne. '')

      ! Option is just a single word (followed by another option or end)
      if( prec_opts(i)(1:1) .eq. '-' .and. prec_opts(i+1)(1:1) .eq. '-' .or. &
          prec_opts(i)(1:1) .eq. '-' .and. prec_opts(i+1)(1:1) .eq. '') then
        !debug: print *, 'A:', trim(prec_opts(i))
        call C_Petsc_Options_Set_Value(trim(prec_opts(i)) // C_NULL_CHAR,  &
                                       C_NULL_CHAR)
        i = i + 1

      ! Option is followed by a switch
      else
        !debug: print *, 'B:', trim(prec_opts(i)), ' ', trim(prec_opts(i+1))
        call C_Petsc_Options_Set_Value(trim(prec_opts(i))   // C_NULL_CHAR,  &
                                       trim(prec_opts(i+1)) // C_NULL_CHAR)
        i = i + 2

      end if
    end do
  end if

  !---------------------------!
  !   Set solver tolerances   !
  !---------------------------!
  Pet % miter = miter
  call C_Petsc_Ksp_Set_Tolerances(Pet % ksp,     &
                                  tol,           &  ! PetscReal rtol
                                  tol,           &  ! PetscReal abstol
                                  Pet % miter)
  !-----------!
  !   Solve   !
  !- - - - - -+------------------------------------------------------!
  !   (This step takes roughly half of the time of filling up the    !
  !    PETSc matrix and setiting the preconditioner, put together)   !
  !------------------------------------------------------------------!
  call C_Petsc_Ksp_Solve(Pet % ksp, Pet % b, Pet % x)

  ! Check if converged
  call C_Petsc_Ksp_Converged_Reason(Pet % ksp, Pet % reason)

  ! Fetch the performed number of iterations
  call C_Petsc_Ksp_Get_Iteration_Number(Pet % ksp, Pet % niter)
  niter = Pet % niter

  ! Fetch the performed number of iterations
  call C_Petsc_Ksp_Get_Residual_Norm(Pet % ksp, fin_res)

  !-----------------------------------------------!
  !   Copy the solution back to T-Flows' vector   ! (A % glo starts from zero)
  !-----------------------------------------------!
  if(Pet % reason > 0 .or.  &            ! converged
     Pet % reason .eq. OUT_OF_ITS) then  ! simply ran out of iterations
    do i = 1, Pet % m_lower
      call C_Petsc_Vec_Get_Values(Pet % x,     &
                                  1,           &
                                  A % glo(i),  &
                                  x(i))
    end do
  else
    ! if(First_Proc()) print *, ' # Warning: linear system failed to converge!'
  end if

  end subroutine

