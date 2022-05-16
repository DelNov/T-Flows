!==============================================================================!
  subroutine Swarm_Mod_Grad_Modeled_Flow(swarm, k)
!------------------------------------------------------------------------------!
!   Stores gradients of modeled Flow parameters for swarm SGS models           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: grid
  type(Turb_Type),  pointer :: Turb
  integer                   :: c
!==============================================================================!

  ! Take aliases for Flow
  Flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  Turb => swarm % pnt_turb

  ! Gradients of turbulent quantities
  call Flow % Grad_Variable(Turb % kin)
  call Flow % Grad_Variable(Turb % eps)

  ! Array for v^2
  do c = -grid % n_bnd_cells, grid % n_cells
    swarm % v2_mod(c) = sqrt(Turb % zeta % n(c) * Turb % kin % n(c))
  end do

  ! Storing the gradients of v2_mod in Work_Mod arrays
  ! (dv^2/dx, dv^2/dy, dv^2/dz)
  call Flow % Grad_Component(swarm % v2_mod, 1, swarm % v2_mod_x)
  call Flow % Grad_Component(swarm % v2_mod, 2, swarm % v2_mod_y)
  call Flow % Grad_Component(swarm % v2_mod, 3, swarm % v2_mod_z)

  end subroutine
