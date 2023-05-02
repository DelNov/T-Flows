!==============================================================================!
  subroutine Solve_Petsc(Pet,                      &
                         solver, prec, prec_opts,  &
                         A, x, b,                  &
                         miter, niter,             &
                         tol, fin_res,             &
                         norm)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)    :: Pet
  character(*)         :: solver          ! solver
  character(*)         :: prec            ! preconditioner
  character(SL)        :: prec_opts(MSI)  ! preconditioner options
  type(Matrix_Type)    :: A
  real                 :: x(-Pet % pnt_grid % n_bnd_cells :  &
                             Pet % pnt_grid % n_cells)
  real                 :: b( Pet % pnt_grid % n_cells)
  integer, intent(in)  :: miter
  integer, intent(out) :: niter
  real,    intent(in)  :: tol      ! tolerance
  real,    intent(out) :: fin_res  ! final residual
  real,    optional    :: norm     ! normalization
!-----------------------------------[Locals]-----------------------------------!
  integer        :: i, j, k
  character(SL)  :: solvers  ! fortran string to store solver
  character(SL)  :: precs    ! fortran string to store preconditioner
  integer        :: l        ! length of the two strings above
!==============================================================================!

  !-----------------------------------------------------------!
  !   Fill up PETSc matrix with values from original matrix   !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !   (This can be very slow, comparable to call to solver)   !
  !-----------------------------------------------------------!
  do i = 1, Pet % m_lower
    do j = A % row(i), A % row(i+1)-1
      k = A % col(j)
      call C_Petsc_Mat_Set_Value(Pet % A,           &  ! matrix
                                 Pet % glo(i),      &  ! row
                                 Pet % glo(k),      &  ! column
                                 A % val(j))           ! matrix entry
    end do
  end do

  ! The following two calls are needed after the calls to MatSetValue
  call C_Petsc_Assemble_Mat(Pet % A)

  !---------------------!
  !   Fill up vectors   ! (Pet % glo starts from zero)
  !---------------------!
  do i = 1, Pet % m_lower
    call C_Petsc_Vec_Set_Value(Pet % x, Pet % glo(i), x(i))
    call C_Petsc_Vec_Set_Value(Pet % b, Pet % glo(i), b(i))
  end do

  ! The following two calls are needed after the calls to VecSetValue
  call C_Petsc_Assemble_Vec(Pet % x)
  call C_Petsc_Assemble_Vec(Pet % b)

  !-----------------------------------!
  !   Set solver and preconditioner   !
  !- - - - - - - - - - - - - - - - - -!
  !   (This can be very slow, much    !
  !    slower than call to solver.)   !
  !-----------------------------------!
  solvers = solver;   l = len_trim(solvers);  solvers(l+1:l+1) = c_null_char
  precs   = prec;     l = len_trim(precs);    precs  (l+1:l+1) = c_null_char
  call C_Petsc_Set_Solver_And_Preconditioner(Pet % ksp,  &  ! solver
                                             Pet % pc,   &  ! preconditioner
                                             Pet % A,    &
                                             solvers,    &
                                             precs)

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
        call C_Petsc_Options_Value(trim(prec_opts(i)), "")
        i = i + 1

      ! Option is followed by a switch
      else
        !debug: print *, 'B:', trim(prec_opts(i)), ' ', trim(prec_opts(i+1))
        call C_Petsc_Options_Value(trim(prec_opts(i)), trim(prec_opts(i+1)))
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
  !   Copy the solution back to T-Flows' vector   ! (Pet % glo starts from zero)
  !-----------------------------------------------!
  if(Pet % reason > 0 .or.  &            ! converged
     Pet % reason .eq. OUT_OF_ITS) then  ! simply ran out of iterations
    do i = 1, Pet % m_lower
      call C_Petsc_Vec_Get_Values(Pet % x,       &
                                  1,             &
                                  Pet % glo(i),  &
                                  x(i))
    end do
  else
    ! if(First_Proc()) print *, ' # Warning: linear system failed to converge!'
  end if

  end subroutine

