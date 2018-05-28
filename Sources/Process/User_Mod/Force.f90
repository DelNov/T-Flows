!==============================================================================!
  subroutine User_Mod_Force(grid, ui, a_matrix, b_vector)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for velocity.      !
!   It is called from "Compute_Velocity" function, just before calling the     !
!   linear solver.  Both system matrix ("a_matrix") and right hand side        !
!   vector ("b_vector") are sent.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Var_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)    :: grid
  type(Var_Type)     :: ui        ! velocity component
  type(Matrix_Type)  :: a_matrix  ! system matrix
  real, dimension(:) :: b_vector  ! right hand side vector
!------------------------------------------------------------------------------!

  !----------------------------------------------------! 
  !                                                    !
  !   Set source depending on the velocity component   !
  !                                                    !
  !----------------------------------------------------! 

  !-------------------------------------------------------!
  !  Set source for velocity component in "x" direction   !
  !-------------------------------------------------------!
  if( ui % name .eq. 'U' ) then  

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

  end subroutine  ! fourth level comments
