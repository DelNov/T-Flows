!==============================================================================!
  subroutine User_Mod_Force(Flow, Por, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This procedure intrduce centrifuge-like force into the T-junction to       !
!   prevent the flow going out through one leg only.  It seemed to have        !
!   helped when it was first introduced on August 21, 2022.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)    :: Flow
  type(Porosity_Type) :: Por
  type(Var_Type)      :: ui        ! velocity component
  type(Matrix_Type)   :: a_matrix  ! system matrix
  real, dimension(:)  :: b_vector  ! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  !----------------------------------------------------!
  !                                                    !
  !   Set source depending on the velocity component   !
  !                                                    !
  !----------------------------------------------------!

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'U' ) then
    do c = Cells_In_Domain_And_Buffers()
      b_vector(c) = b_vector(c) + 1.0e4 * Grid % xc(c) * Grid % vol(c)
    end do
  end if

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'V' ) then

  end if

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'W' ) then

  end if

  end subroutine
