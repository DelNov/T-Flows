!==============================================================================!
  subroutine Front_Mod_Place_At_Var_Value(front, phi, sol, phi_e, verbose)
!------------------------------------------------------------------------------!
!   Places surface where variable phi has value phi_e                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n   => r_node_01,  &  ! value at the static grid nodes
                      phi_c   => r_cell_01,  &  ! cell values of phi
                      phi_o   => r_cell_02,  &  ! original values of phi
                      phi_cen => r_cell_03,  &
                      phi_src => r_cell_04
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type),   target :: front
  type(Var_Type),    target :: phi
  type(Solver_Type), target :: sol   ! needed for smoothing
  real                      :: phi_e
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Field_Type),  pointer :: flow
  type(Vert_Type),   pointer :: vert(:)
  type(Elem_Type),   pointer :: elem(:)
  type(Matrix_Type), pointer :: a
  integer,           pointer :: nv, ne
  integer, allocatable       :: n_cells_v(:)
  integer                    :: c, c1, c2, s, j, n1, n2, run, nb, nc, n, nn
  integer                    :: v, n_vert, n_verts_in_buffers
  integer                    :: en(12,2)  ! edge numbering
  real                       :: phi1, phi2, xn1, yn1, zn1, xn2, yn2, zn2, w1, w2
  real                       :: surf_v(3), dist
!------------------------------------------------------------------------------!
  include 'Front_Mod/Edge_Numbering.f90'
!==============================================================================!

  ! Take aliases
  grid => front % pnt_grid
  flow => front % pnt_flow
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem
  a    => sol % a
  nb   = grid % n_bnd_cells
  nc   = grid % n_cells
  nn   = grid % n_nodes

  call Front_Mod_Initialize(front)
  call Field_Mod_Interpolate_Cells_To_Nodes(flow, phi % n, phi_n(1:nn))

  CALL GRID_MOD_SAVE_DEBUG_VTU(GRID, 'PHI_N',       &
                               SCALAR_NODE=PHI_N,   &
                               SCALAR_NAME='PHI_N')

  !-----------------------------------------------!
  !   This is a bit ad-hoc - if some points are   !
  !   "exactly" 0.5, increase them a little bit   !
  !-----------------------------------------------!
  do n = 1, nn
    if(Math_Mod_Approx_Real(phi_n(n), 0.5)) phi_n(n) = 0.5 + MILI
  end do

  !-----------------------------!
  !   Smooth the VOF function   !
  !-----------------------------!

  ! Copy the values from phi % n to local variable
  phi_o(-nb:nc) = phi % n(-nb:nc)
  phi_c(-nb:nc) = phi % n(-nb:nc)

  do run = 1, 16
    phi_src(:) = 0.0
    phi_cen(:) = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 > 0) then
        phi_src(c1) = phi_src(c1) + A % fc(s) * phi_c(c2)
        phi_cen(c1) = phi_cen(c1) + A % fc(s)
        phi_src(c2) = phi_src(c2) + A % fc(s) * phi_c(c1)
        phi_cen(c2) = phi_cen(c2) + A % fc(s)
      end if
    end do
    do c = 1, grid % n_cells
      phi_c(c) = 0.1 * phi_c(c) + 0.9 * phi_src(c) / phi_cen(c)
    end do
  end do

  phi % n(:) = phi_c(:)
  call Field_Mod_Grad_Variable(flow, phi)

  allocate(n_cells_v(grid % n_cells))
  n_cells_v(:) = 0

  ! Find vertices in each cell
  nv = 0
  ne = 0

  !------------------------------------------------------------!
  !                                                            !
  !   Browse through all cells in search of surface vertices   !
  !                                                            !
  !------------------------------------------------------------!
  do c = 1, grid % n_cells

    n_vert = 0

    ! Fetch the edges for this cell
    if(grid % cells_n_nodes(c) .eq. 4) en = edg_tet
    if(grid % cells_n_nodes(c) .eq. 5) en = edg_pyr
    if(grid % cells_n_nodes(c) .eq. 6) en = edg_wed
    if(grid % cells_n_nodes(c) .eq. 8) en = edg_hex

    !------------------------------------------------------!
    !   Browse through edges to find intersection points   !
    !------------------------------------------------------!
    do j = 1, 12  ! max number of edges
      n1 = grid % cells_n( en(j,1), c )
      n2 = grid % cells_n( en(j,2), c )

      phi1 = phi_n(n1)
      phi2 = phi_n(n2)

      ! There is a vertex between these two edges
      if( ((phi2-phi_e) * (phi_e-phi1)) >= MICRO ) then
        n_cells_v(c) = n_cells_v(c) + 1

        nv = nv + 1

        n_vert = n_vert + 1

        xn1 = grid % xn(n1)
        yn1 = grid % yn(n1)
        zn1 = grid % zn(n1)

        xn2 = grid % xn(n2)
        yn2 = grid % yn(n2)
        zn2 = grid % zn(n2)

        w1 = abs(phi2-phi_e) / abs(phi2-phi1)
        w2 = 1.0 - w1

        ! All vertices have to be stored
        vert(nv) % x_n = xn1*w1 + xn2*w2
        vert(nv) % y_n = yn1*w1 + yn2*w2
        vert(nv) % z_n = zn1*w1 + zn2*w2

      end if

    end do  ! through edges

    !---------------------------------!
    !   Some points have been found   !
    !---------------------------------!
    if(n_vert > 2) then

      ! Surface vector
      surf_v(1) = phi % x(c)
      surf_v(2) = phi % y(c)
      surf_v(3) = phi % z(c)

      ! If valid elements were formed
      if(n_vert .eq. 3) call Front_Mod_Handle_3_Points(front, surf_v)
      if(n_vert .eq. 4) call Front_Mod_Handle_4_Points(front, surf_v)
      if(n_vert .eq. 5) call Front_Mod_Handle_5_Points(front, surf_v)
      if(n_vert .eq. 6) call Front_Mod_Handle_6_Points(front, surf_v)
      if(n_vert .eq. 7) then
        print *, '# ERROR: seven vertices in an intersection!'
        stop

      end if

      ! Store at which cell the surface resides
      elem(ne) % cell = c

    end if

  end do
  if(verbose) then
    print *, '# Cummulative number of elements found: ', ne
    print *, '# Cummulative number of vertices found: ', nv
  end if

  !-----------------------!
  !                       !
  !   Compress vertices   !
  !                       !
  !-----------------------!
  call Front_Mod_Compress_Vertices(front, verbose)

  !----------------!
  !                !
  !   Find sides   !
  !                !
  !----------------!
  call Front_Mod_Find_Connectivity(front, verbose)

  !--------------------------------!
  !   Find nearest cell and node   !
  !--------------------------------!
  n_verts_in_buffers = 0
  do v = 1, nv
    call Front_Mod_Find_Nearest_Cell(front, v, n_verts_in_buffers)
    call Front_Mod_Find_Nearest_Node(front, v)
  end do

  !--------------------------------------!
  !                                      !
  !   Calculate geometrical quantities   !
  !                                      !
  !--------------------------------------!

  ! Element geometry has changed, recompute geometrical quantities
  call Front_Mod_Find_Vertex_Elements(front)
  call Front_Mod_Calculate_Element_Centroids(front)
  call Front_Mod_Calculate_Element_Normals(front, phi)

  ! Restore the true values of phi
  phi % n(:) = phi_o(:)

  return

  end subroutine
