!==============================================================================!
  subroutine User_Mod_Source(Grid, Flow, Turb, phi, sc)
!------------------------------------------------------------------------------!
!   User-defined source terms for scalar variables.                            !
!   Be mindful that this can be performed on the GPUs                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),     target   :: Grid
  type(Field_Type),    target   :: Flow
  type(Turb_Type),     target   :: Turb
  type(Var_Type),      target   :: phi   !! scalar being solved
  integer, intent(in), optional :: sc    !! scalar rank (optional)
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:)  ! pointer to matrix values
  integer,   contiguous, pointer :: dia(:)  ! pointer to matrix diagonal
  real,      contiguous, pointer :: b(:)    ! pointer to right hand side
  integer                        :: c       ! cell counter
!==============================================================================!

  !-----------------------------!
  !   First take some aliases   !
  !-----------------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % A % dia
  b   => Flow % Nat % b

  !------------------------------------!
  !                                    !
  !   Insert some source for scalars   !
  !                                    !
  !------------------------------------!
  if(present(sc)) then

    !-----------------------!
    !   Focus on scalar 1   !
    !-----------------------!
    if(sc .eq. 1) then
      !$tf-acc loop begin
      do c = Cells_In_Domain()
        b(c) = b(c) + 0.0 * Grid % vol(c)
      end do
      !$tf-acc loop end
    end if

  end if

  end subroutine
