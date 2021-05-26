!==============================================================================!
  subroutine Improve_Mesh_Quality(Surf, smooth, phi_e, verbose)
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
  real                      :: phi_e
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Vert_Type),   pointer :: Vert(:)
  type(Elem_Type),   pointer :: Elem(:)
  integer,           pointer :: nv, ne
  integer, allocatable       :: n_cells_v(:)
  integer                    :: c, c1, c2, s, j, n1, n2, run, nb, nc, nn
  integer                    :: v, n_vert, n_verts_in_buffers
  integer                    :: en(12,2)  ! edge numbering
  real                       :: phi1, phi2, xn1, yn1, zn1, xn2, yn2, zn2, w1, w2
  real                       :: surf_v(3), dist
!------------------------------------------------------------------------------!
  include 'Surf_Mod/Edge_Numbering.f90'
!==============================================================================!

  ! Take aliases
  Grid => Surf % pnt_grid
  Flow => Surf % pnt_flow
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  nn   =  Grid % n_nodes

  !--------------------------------------!
  !   Relax topology, followed by mesh   !
  !    smoothing, in a few iterations    !
  !--------------------------------------!
  do j = 1, 3
    call Surf % Relax_Topology()
    call Surf % Smooth(smooth, phi_e)
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

  ! The rest is still experimental
  call Surf % Refine(4)
  do j = 1, 3
    call Surf % Relax_Topology()
    call Surf % Smooth(smooth, phi_e)
  end do
  do j = 1, 3
    call Surf % Relax_Geometry()
    call Surf % Smooth(smooth, phi_e)
  end do

  end subroutine
