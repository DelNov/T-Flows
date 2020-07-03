!==============================================================================!
  subroutine Swarm_Mod_Grad_Modeled_Flow(swarm, turb, k)
!------------------------------------------------------------------------------!
!   Stores gradients of modeled flow parameters for swarm SGS models           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: w_mod_x => r_cell_01,   &
                      w_mod_y => r_cell_02,   &
                      w_mod_z => r_cell_03,   &
                      w_mod   => r_cell_04
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
  integer                  :: k
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  integer                      :: c, nb, nc
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  nb = grid % n_bnd_cells
  nc = grid % n_cells

  ! Gradients of turbulent quantities
  call Field_Mod_Grad_Variable(flow, turb % kin)
  call Field_Mod_Grad_Variable(flow, turb % eps)

  ! Array for v^2
  do c = -grid % n_bnd_cells, grid % n_cells
    w_mod(c) = sqrt(turb % zeta % n(c) * turb % kin % n(c))
  end do

  ! Storing the gradients of w_mod in Work_Mod arrays
  ! (dv^2/dx, dv^2/dy, dv^2/dz)
  call Field_Mod_Grad_Component(flow, w_mod(-nb:nc), 1, w_mod_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, w_mod(-nb:nc), 2, w_mod_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, w_mod(-nb:nc), 3, w_mod_z(-nb:nc))

  ! Wall-normal velocity contribution from the modeled part (zeta)
  do c = -grid % n_bnd_cells, grid % n_cells
    swarm % w_mod(c) = w_mod(c)
    swarm % w_x(c)   = w_mod_x(c)
    swarm % w_y(c)   = w_mod_y(c)
    swarm % w_z(c)   = w_mod_z(c)
  end do

  end subroutine
