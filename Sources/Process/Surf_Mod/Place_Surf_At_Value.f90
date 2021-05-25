!==============================================================================!
  subroutine Place_Surf_At_Value(Surf, sharp, smooth, phi_e, verbose)
!------------------------------------------------------------------------------!
!   Places surface where variable phi has value phi_e                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n => r_node_01  ! value at the static Grid nodes
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf
  type(Var_Type),    target :: sharp
  type(Var_Type),    target :: smooth
  real                      :: phi_e
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Vert_Type),   pointer :: Vert(:)
  type(Elem_Type),   pointer :: Elem(:)
  type(Matrix_Type), pointer :: A
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

  call Cpu_Timer % Start('Creating_Surface_From_Vof_Function')

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

  call Surf % Initialize_Surf()
  call Flow % Interpolate_Cells_To_Nodes(sharp % n, phi_n(1:nn))

  !---------------------------------------------!
  !   Take gradients from the smooth function   !
  !   These might be already computed - check   !
  !---------------------------------------------!
  call Flow % Grad_Variable(smooth)

  allocate(n_cells_v(Grid % n_cells))
  n_cells_v(:) = 0

  ! Find vertices in each cell
  nv = 0
  ne = 0

  !------------------------------------------------------------!
  !                                                            !
  !   Browse through all cells in search of surface vertices   !
  !                                                            !
  !------------------------------------------------------------!
  do c = 1, Grid % n_cells

    n_vert = 0

    ! Fetch the edges for this cell
    if(Grid % cells_n_nodes(c) .eq. 4) en = edg_tet
    if(Grid % cells_n_nodes(c) .eq. 5) en = edg_pyr
    if(Grid % cells_n_nodes(c) .eq. 6) en = edg_wed
    if(Grid % cells_n_nodes(c) .eq. 8) en = edg_hex

    !------------------------------------------------------!
    !   Browse through edges to find intersection points   !
    !------------------------------------------------------!
    do j = 1, 12  ! max number of edges
      n1 = Grid % cells_n( en(j,1), c )
      n2 = Grid % cells_n( en(j,2), c )

      phi1 = phi_n(n1)
      phi2 = phi_n(n2)

      ! There is a vertex between these two edges
      if( ((phi2-phi_e) * (phi_e-phi1)) >= MICRO ) then
        n_cells_v(c) = n_cells_v(c) + 1

        nv = nv + 1

        n_vert = n_vert + 1

        xn1 = Grid % xn(n1)
        yn1 = Grid % yn(n1)
        zn1 = Grid % zn(n1)

        xn2 = Grid % xn(n2)
        yn2 = Grid % yn(n2)
        zn2 = Grid % zn(n2)

        w1 = abs(phi2-phi_e) / abs(phi2-phi1)
        w2 = 1.0 - w1

        ! All vertices have to be stored
        Vert(nv) % x_n = xn1*w1 + xn2*w2
        Vert(nv) % y_n = yn1*w1 + yn2*w2
        Vert(nv) % z_n = zn1*w1 + zn2*w2

      end if

    end do  ! through edges

    !---------------------------------!
    !   Some points have been found   !
    !---------------------------------!
    if(n_vert > 0) then

      ! Surface vector
      surf_v(1) = smooth % x(c)
      surf_v(2) = smooth % y(c)
      surf_v(3) = smooth % z(c)

      ! If valid elements were formed (last argument: enforce_triangles)
      if(n_vert .eq. 3) call Surf % Handle_3_Points(surf_v, .true.)
      if(n_vert .eq. 4) call Surf % Handle_4_Points(surf_v, .true.)
      if(n_vert .eq. 5) call Surf % Handle_5_Points(surf_v, .true.)
      if(n_vert .eq. 6) call Surf % Handle_6_Points(surf_v, .true.)
      if(n_vert .eq. 7) then
        print *, '# ERROR: seven vertices in an intersection!'
        stop
      end if
    end if

  end do
  if(verbose) then
    print '(a40,i8)', ' # Cummulative number of elements found:', ne
    print '(a40,i8)', ' # Cummulative number of vertices found:', nv
  end if

  !-----------------------!
  !                       !
  !   Compress vertices   !
  !                       !
  !-----------------------!
  call Surf % Compress_Vertices(verbose)

  !---------------------------------!
  !                                 !
  !   Find and check connectivity   !
  !                                 !
  !---------------------------------!
  call Surf % Find_Sides(verbose)          ! Front calls the same here
  call Surf % Find_Surf_Elements(verbose)  ! Front calls Find_Front_Elements
  call Surf % Check_Elements(verbose)      ! Front calls the same here

  !--------------------------------!
  !   Find nearest cell and node   !
  !--------------------------------!
  n_verts_in_buffers = 0
  do v = 1, nv
    call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
    call Surf % Vert(v) % Find_Nearest_Node()
  end do

  !-------------------------------!
  !                               !
  !   Calculate element normals   !
  !                               !
  !-------------------------------!
  call Surf % Calculate_Element_Centroids()
  call Surf % Calculate_Element_Normals(smooth)

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


  call Cpu_Timer % Stop('Creating_Surface_From_Vof_Function')

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
