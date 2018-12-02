!==============================================================================!
  subroutine Field_Mod_Allocate(flow, grid)
!------------------------------------------------------------------------------!
!   Allocates memory for variables. It is called either from LoaRes            !
!   or from Processor.                                                         !
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)         :: flow
  type(Grid_Type),  target :: grid
!==============================================================================!

  ! Store the pointer to a grid
  flow % pnt_grid => grid

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Allocate_Solution('U', flow % u, grid)
  call Var_Mod_Allocate_Solution('V', flow % v, grid)
  call Var_Mod_Allocate_Solution('W', flow % w, grid)

  ! ... and their gradients
  call Var_Mod_Allocate_Gradients(flow % u)
  call Var_Mod_Allocate_Gradients(flow % v)
  call Var_Mod_Allocate_Gradients(flow % w)

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only('PP', flow % pp, grid)
  call Var_Mod_Allocate_New_Only('P',  flow % p,  grid)

  ! Pressure gradients are needed too
  call Var_Mod_Allocate_Gradients(flow % p)

  ! Mass flow rates at cell faces are always needed
  allocate(flow % flux(grid % n_faces));  flow % flux = 0.

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(heat_transfer) then
    call Var_Mod_Allocate_Solution('T', flow % t, grid)
  end if ! heat_transfer

  end subroutine
