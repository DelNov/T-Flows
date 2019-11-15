!==============================================================================!
  subroutine User_Mod_Source(flow, phi, a, b)
!------------------------------------------------------------------------------!
!   This is a prototype of a function for customized source for scalar.        !
!   It is called from "Compute_Scalar" function, just before calling the       !
!   linear solver.  Both system matrix ("a") and right hand side        !
!   vector ("b") are sent should the user want to stabilize the         !
!   system for always positive variables, for example.                         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,   only: Grid_Type
  use Var_Mod,    only: Var_Type
  use Field_Mod
  use Matrix_Mod, only: Matrix_Type
  use Bulk_Mod,   only: Bulk_Type
  use Const_Mod,  only: PI
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Var_Type)           :: phi
  type(Matrix_Type)        :: a
  real, dimension(:)       :: b
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t
  type(Bulk_Type), pointer :: bulk
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  bulk => flow % bulk

  !-----------------------------------------------------! 
  !                                                     !
  !   Set source depending on which variable you have   !
  !                                                     !
  !-----------------------------------------------------! 

  !-------------------------------!
  !  Set source for temperature   !
  !-------------------------------!
  if( phi % name .eq. 'T' ) then

    call Comm_Mod_Global_Sum_Real(heat)
    call Comm_Mod_Global_Sum_Real(heated_area)
    heat_flux = heat / (heated_area + TINY)

    do c = 1, grid % n_cells
      b(c) = b(c) - 2.0 * pi * heat_flux * w % n(c)  &
           / bulk % flux_z * grid % vol(c)
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
