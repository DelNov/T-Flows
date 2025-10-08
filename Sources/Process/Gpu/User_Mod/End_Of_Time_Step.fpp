!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Grid, Flow, Turb)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),     target   :: Grid
  type(Field_Type),    target   :: Flow
  type(Turb_Type),     target   :: Turb
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: b(:)    ! pointer to right hand side
  integer                        :: c       ! cell counter
!==============================================================================!

  !-----------------------------!
  !   First take some aliases   !
  !-----------------------------!
  b => Flow % Nat % b

  !-------------------------------------!
  !   Pretend you do something useful   !
  !-------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain()
    b(c) = b(c) + 0.0 * Grid % vol(c)
  end do
  !$tf-acc loop end

  end subroutine
