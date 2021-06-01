!==============================================================================!
  subroutine Improve_Mesh_Quality(Surf, smooth, verbose)
!------------------------------------------------------------------------------!
!   Improves the mesh quality for the surface                                  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n => r_node_01  ! value at the static Grid nodes
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf
  type(Var_Type),    target :: smooth
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  integer :: j
!==============================================================================!

  !--------------------------------------!
  !   Relax topology, followed by mesh   !
  !    smoothing, in a few iterations    !
  !--------------------------------------!
  do j = 1, 3
    call Surf % Relax_Topology()
    call Surf % Smooth_Surf(smooth)
  end do

  !-> ! From this point ...
  !-> do v = 1, nv
  !->   dist = norm2( (/Surf % Vert(v) % x_n,  &
  !->                   Surf % Vert(v) % y_n,  &
  !->                   Surf % Vert(v) % z_n/) )
  !->   Surf % Vert(v) % x_n = Surf % Vert(v) % x_n * 0.25 / dist
  !->   Surf % Vert(v) % y_n = Surf % Vert(v) % y_n * 0.25 / dist
  !->   Surf % Vert(v) % z_n = Surf % Vert(v) % z_n * 0.25 / dist
  !-> end do
  !-> n_verts_in_buffers = 0
  !-> do v = 1, nv
  !->   call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
  !->   call Surf % Vert(v) % Find_Nearest_Node()
  !-> end do
  !-> ! ... down to here is just for development

  ! Element geometry has changed, recompute geometrical quantities
  call Surf % Find_Vertex_Elements()
  call Surf % Calculate_Element_Centroids()
  call Surf % Calculate_Element_Normals(smooth)

  return

  !-> ! The rest is still experimental
  !-> call Surf % Refine(4)
  !-> do j = 1, 3
  !->   call Surf % Relax_Topology()
  !->   call Surf % Smooth(smooth)
  !-> end do
  !-> do j = 1, 3
  !->   call Surf % Relax_Geometry()
  !->   call Surf % Smooth(smooth)
  !-> end do

  end subroutine
