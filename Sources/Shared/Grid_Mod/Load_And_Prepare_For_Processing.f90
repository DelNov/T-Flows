!==============================================================================!
  subroutine Load_And_Prepare_For_Processing(Grid, d)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid
  integer,          intent(in)    :: d
!==============================================================================!

  call Control % Read_Problem_Name(problem_name(d))

  ! Load the finite volume Grid
  call Grid % Load_Cfn((/This_Proc(), N_Procs()/), domain=d)
  call Grid % Load_Dim((/This_Proc(), N_Procs()/), domain=d)

  ! Determine threads for OpenMP runs
  call Control % Max_Threads(Grid % Vect % d_threads, .true.)
  call Grid % Determine_Threads()

  call Grid % Calculate_Face_Geometry()

  ! Find communication patterns
  call Grid % Form_Cells_Comm()
  call Grid % Form_Maps_For_Backup()

  end subroutine

