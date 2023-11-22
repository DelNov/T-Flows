!==============================================================================!
  subroutine Grad_Modeled_Flow(Swarm)
!------------------------------------------------------------------------------!
!   Stores gradients of modeled Flow parameters for swarm SGS models           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Turb_Type),  pointer :: Turb
  integer                   :: c
!==============================================================================!

  ! Take aliases for Flow
  Flow => Swarm % pnt_flow
  Grid => Swarm % pnt_grid
  Turb => Swarm % pnt_turb

  ! Gradients of turbulent quantities
  call Flow % Grad_Variable(Turb % kin)
  call Flow % Grad_Variable(Turb % eps)

  ! Array for v^2
  do c = -Grid % n_bnd_cells, Grid % n_cells
    Swarm % v2_mod(c) = sqrt(Turb % zeta % n(c) * Turb % kin % n(c))
  end do

  ! Storing the gradients of v2_mod in Work_Mod arrays
  ! (dv^2/dx, dv^2/dy, dv^2/dz)
  call Flow % Grad_Component(Swarm % v2_mod, 1, Swarm % v2_mod_x)
  call Flow % Grad_Component(Swarm % v2_mod, 2, Swarm % v2_mod_y)
  call Flow % Grad_Component(Swarm % v2_mod, 3, Swarm % v2_mod_z)

  end subroutine
