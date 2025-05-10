!==============================================================================!
  subroutine User_Mod_Source(Flow, phi, a, b)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a") and right hand side               !
!   vector ("b") are sent should the user want to stabilize the                !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Var_Type),   target :: phi
  type(Matrix_Type)        :: a
  real, dimension(:)       :: b
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  integer                  :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name .eq. 'T' ) then
    do c = Cells_In_Domain_And_Buffers()
      b(c) = b(c) - w % n(c) / bulk % w          &
                  * Flow % heat / bulk % area_z  &
                  * Grid % vol(c)
    end do
  end if

  !---------------------------------------------!
  !  Set source for turbulent kintetic energy   !
  !---------------------------------------------!
  if( phi % name .eq. 'KIN' ) then

  end if

  !---------------------------!
  !  Set source for epsilon   !
  !---------------------------!
  if( phi % name .eq. 'EPS' ) then

  end if

  end subroutine  ! fourth level comments
