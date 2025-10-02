!==============================================================================!
  subroutine Improve_Mesh_Quality(Surf, smooth)
!------------------------------------------------------------------------------!
!>  The Improve_Mesh_Quality subroutine focuses on optimizing the mesh
!>  structure of a surface object. This optimization is achieved through a
!>  series of relaxation and smoothing operations that improve the topological
!>  and geometrical characteristics of the mesh.  This strategy of improving
!>  the mesh quality is taken from TRIPOS (https://github.com/Niceno/TRIPOS)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf    !! parent class
  type(Var_Type),   target :: smooth  !! smooth variant of the scalar field
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
  call Surf % Calculate_Element_Normals()

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
