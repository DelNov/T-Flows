!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         blend_matrix)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)          :: Pet
  character(*),  intent(in)  :: solver          ! solver
  character(*),  intent(in)  :: prec            ! preconditioner
  character(SL), intent(in)  :: prec_opts(MSI)  ! preconditioner options
  type(Matrix_Type)          :: A
  real                       :: x(-Pet % pnt_grid % n_bnd_cells :  &
                                   Pet % pnt_grid % n_cells)
  real                       :: b( Pet % pnt_grid % n_cells)
  integer,       intent(in)  :: miter
  integer,       intent(out) :: niter
  real,          intent(in)  :: tol      ! tolerance
  real,          intent(out) :: fin_res  ! final residual
  logical,       intent(in)  :: blend_matrix
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i, j, k
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
        call C_Petsc_Mat_Set_Value(Pet % A,     &  ! matrix
                                   A % glo(i),  &  ! row
                                   A % glo(k),  &  ! column
                                   A % val(j))     ! matrix entry
      end do
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
    do while(i < MSI .and. prec_opts(i)(1:1) .ne. '')

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

