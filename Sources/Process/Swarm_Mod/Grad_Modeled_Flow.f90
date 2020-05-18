!==============================================================================!
  subroutine Swarm_Mod_Grad_Modeled_Flow(swarm, turb, k)
!------------------------------------------------------------------------------!
!   Stores gradients of modeled flow parameters for swarm SGS models           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: kin_x         => r_cell_01,   &
                      kin_y         => r_cell_02,   &
                      kin_z         => r_cell_03,   &
                      zeta_x        => r_cell_04,   &
                      zeta_y        => r_cell_05,   &
                      zeta_z        => r_cell_06,   &
                      kin           => r_cell_07,   &
                      zeta          => r_cell_08
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  type(Turb_Type),  target :: turb
  integer                  :: k
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  integer                   :: c, c2, cc               ! nearest cell
  integer                   :: nb, nc                  ! nearest cell
  real                      :: up, vp, wp              ! velocity at particle
  real                      :: w2_dx2, w2_dy2, w2_dz2  ! v2 quantity components
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid
  nc   =  grid % n_cells
  nb   =  grid % n_bnd_cells

  ! Arrays for modeled quantities (needed if Hybrid_Les_Rans model is used)
  do c = -grid % n_bnd_cells, grid % n_cells
    kin(c)  = turb % kin_mean(c)
    zeta(c) = turb % zeta_mean(c)
  end do

  ! Storing the gradients of kin_mean and zeta_mean in Work_Mod arrays
  call Field_Mod_Grad_Component(flow, kin (-nb:nc), 1, kin_x (-nb:nc))
  call Field_Mod_Grad_Component(flow, kin (-nb:nc), 2, kin_y (-nb:nc))
  call Field_Mod_Grad_Component(flow, kin (-nb:nc), 3, kin_z (-nb:nc))
  call Field_Mod_Grad_Component(flow, zeta(-nb:nc), 1, zeta_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, zeta(-nb:nc), 2, zeta_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, zeta(-nb:nc), 3, zeta_z(-nb:nc))

  do c = -grid % n_bnd_cells, grid % n_cells

    ! v2 quantity components
    w2_dx2 = kin_x(c) * zeta_x(c)
    w2_dy2 = kin_y(c) * zeta_y(c)
    w2_dz2 = kin_z(c) * zeta_z(c)

    ! d(zeta)/dx/dx
    swarm % w_mod_x(c) = abs(w2_dx2)
    swarm % w_mod_y(c) = abs(w2_dy2)
    swarm % w_mod_z(c) = abs(w2_dz2)

    ! Modeled TKE for swarm all over the grid
    swarm % kin_mod(c)  = kin(c)

  end do

  end subroutine
