!==============================================================================!
  subroutine Swarm_Mod_Grad_Modeled_Flow(swarm, turb, k)
!------------------------------------------------------------------------------!
!       Stores gradients of modeled flow parameters for swarm SGS models       !
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
  type(Field_Type),    pointer :: flow
  type(Grid_Type),     pointer :: grid
  integer                      :: c, c2, cc                      ! nearest cell
  real                         :: up, vp, wp                     ! velocity at particle
  real                         :: w2_dx2, w2_dy2, w2_dz2         ! v2 quantity components 
!==============================================================================!

  ! Take aliases for flow
  flow => swarm % pnt_flow
  grid => swarm % pnt_grid

  ! Arrays for modeled quantities (needed if Hybrid_Les_Rans model is used)
  !do c = 1, grid % n_cells
  do c = -grid % n_bnd_cells, grid % n_cells
    kin(c)  = turb % kin_mean(c)
    zeta(c) = turb % zeta_mean(c)
  end do

  ! storing the gradients of kin_mean and zeta_mean in Work_Mod arrays
  call Field_Mod_Grad_Component(flow, kin,  1, kin_x)   ! dk/dx
  call Field_Mod_Grad_Component(flow, kin,  2, kin_y)   ! dk/dy
  call Field_Mod_Grad_Component(flow, kin,  3, kin_z)   ! dk/dz
  call Field_Mod_Grad_Component(flow, zeta, 1, zeta_x)  ! dzeta/dx
  call Field_Mod_Grad_Component(flow, zeta, 2, zeta_y)  ! dzeta/dy
  call Field_Mod_Grad_Component(flow, zeta, 3, zeta_z)  ! dzeta/dz

  !do  c = 1, grid % n_cells
  do c = -grid % n_bnd_cells, grid % n_cells
 
    ! v2 quantity components 
    w2_dx2 = kin_x(c) * zeta_x(c)  
    w2_dy2 = kin_y(c) * zeta_y(c)   
    w2_dz2 = kin_z(c) * zeta_z(c)    
    
    !! not sure if this will be really needed! 
    !sign1 = sign(w2_dx2, w2_dx2) / abs(w2_dx2)
    !sign2 = sign(w2_dy2, w2_dy2) / abs(w2_dy2)
    !sign3 = sign(w2_dz2, w2_dz2) / abs(w2_dz2)

    ! d(zeta)/dx/dx
    swarm % w_mod_x(c) = abs(w2_dx2)
    swarm % w_mod_y(c) = abs(w2_dy2)
    swarm % w_mod_z(c) = abs(w2_dz2)
   
    ! Modeled TKE for swarm all over the grid
    swarm % kin_mod(c)  = kin(c)

  end do 
 
  end subroutine
