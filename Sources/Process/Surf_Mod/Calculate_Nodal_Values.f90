!==============================================================================!
  subroutine Surf_Mod_Calculate_Nodal_Values(surf, phi)
!------------------------------------------------------------------------------!
!   Calculates nodal values of variable phi                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n => r_node_01  ! value at the (static) grid nodes
!------------------------------------------------------------------------------!
!   Be careful with the above variables from Work_Mod.  They are used by       !
!   two subroutines in Surf_Mod, hence values shouldn't be changed elsewhere.  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  type(Var_Type),  target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  logical, allocatable     :: is_node_bnd(:)
  integer                  :: c, ni, n, s
  real                     :: xc, yc, zc, xn, yn, zn
!==============================================================================!

  ! Take aliases
  grid => surf % pnt_grid

  ! Mark boundary nodes (one should set true boundary values there)
  allocate(is_node_bnd(grid % n_nodes))
  is_node_bnd(:) = .false.
  do s = 1, grid % n_faces
    c = grid % faces_c(2, s)
    if(c < 0) then
      do ni = 1, grid % faces_n_nodes(s)  ! 3 or 4
        is_node_bnd( grid % faces_n(ni, s) ) = .true.
      end do
    end if
  end do

  !-------------------------------------------------!
  !   Compute values of variable phi on the nodes   !
  !-------------------------------------------------!

  ! Extrapolate from cells to nodes with variable's gradients
  call Grad_Mod_Variable(phi, boundary = .true.)

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    xc = grid % xc(c)
    yc = grid % yc(c)
    zc = grid % zc(c)
    do ni = 1, grid % cells_n_nodes(c)
      n  = grid % cells_n(ni, c)
      xn = grid % xn(n)
      yn = grid % yn(n)
      zn = grid % zn(n)
      if( .not. is_node_bnd(n) ) then
        phi_n(n) = phi_n(n) + phi % n(c)              &
                            + phi % x(c) * (xn - xc)  &
                            + phi % y(c) * (yn - yc)  &
                            + phi % z(c) * (zn - zc)
      end if
    end do
  end do

  do n = 1, grid % n_nodes
    if( .not. is_node_bnd(n) ) then
      phi_n(n) = phi_n(n) / grid % nodes_n_cells(n)
    end if
  end do

  end subroutine
