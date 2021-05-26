!==============================================================================!
  subroutine Place_Front_At_Value(Front, sharp, smooth, phi_e, verbose)
!------------------------------------------------------------------------------!
!   Places surface where variable phi has value phi_e                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: phi_n => r_node_01  ! value at the static Grid nodes
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  type(Var_Type),    target :: sharp
  type(Var_Type),    target :: smooth
  real                      :: phi_e
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Vert_Type),   pointer :: Vert(:)
  type(Elem_Type),   pointer :: Elem(:)
  integer,           pointer :: nv, ne
  integer                    :: nv_tot, ne_tot
  integer, allocatable       :: n_cells_v(:)
  integer                    :: c, c1, c2, s, j, n1, n2, run, nb, nc, n, nn
  integer                    :: v, n_vert, n_verts_in_buffers, i_nod
  integer                    :: en(12,2)  ! edge numbering
  real                       :: phi1, phi2, xn1, yn1, zn1, xn2, yn2, zn2, w1, w2
  real                       :: surf_v(3), dist
!------------------------------------------------------------------------------!
  include 'Front_Mod/Edge_Numbering.f90'
!==============================================================================!

  call Cpu_Timer % Start('Creating_Front_From_Vof_Function')

  ! Take aliases
  Grid => Front % pnt_grid
  Flow => Front % pnt_flow
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem
  nb   =  Grid % n_bnd_cells
  nc   =  Grid % n_cells
  nn   =  Grid % n_nodes

  call Front % Initialize_Front()
  call Flow % Interpolate_Cells_To_Nodes(sharp % n, phi_n(1:nn))

! call Grid % Save_Debug_Vtu('phi_c',               &
!                            scalar_cell=phi % n,   &
!                            scalar_name='phi_c')
!
! call Grid % Save_Debug_Vtu('phi_n',             &
!                            scalar_node=phi_n,   &
!                            scalar_name='phi_n')

  !-----------------------------------------------!
  !   This is a bit ad-hoc - if some points are   !
  !    "exactly" 0.5, change them a little bit    !
  !-----------------------------------------------!
  do c = 1, nc
    do i_nod = 1, Grid % cells_n_nodes(c)
      n = Grid % cells_n(i_nod, c)
      if(Math % Approx_Real(phi_n(n), 0.5)) then
        if(sharp % n(c) < sharp % o(c)) then
          phi_n(n) = 0.5 - MILI
        else
          phi_n(n) = 0.5 + MILI
        end if
      end if
    end do
  end do

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
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells

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
    if(n_vert > 2) then

      ! Surface vector
      surf_v(1) = smooth % x(c)
      surf_v(2) = smooth % y(c)
      surf_v(3) = smooth % z(c)

      ! If valid elements were formed (last parameter: enforce_triangles)
      if(n_vert .eq. 3) call Front % Handle_3_Points(surf_v, .false.)
      if(n_vert .eq. 4) call Front % Handle_4_Points(surf_v, .false.)
      if(n_vert .eq. 5) call Front % Handle_5_Points(surf_v, .false.)
      if(n_vert .eq. 6) call Front % Handle_6_Points(surf_v, .false.)
      if(n_vert .eq. 7) then
        print *, '# ERROR: seven vertices in an intersection!'
        stop

      end if

      ! Store at which cell the surface resides
      Elem(ne) % cell = c

    end if

  end do

  if(verbose) then
    ne_tot = ne
    nv_tot = nv
    call Comm_Mod_Global_Sum_Int(ne_tot)
    call Comm_Mod_Global_Sum_Int(nv_tot)
    if(this_proc < 2) then
      print '(a40,i8)', ' # Cummulative number of elements found:', ne_tot
      print '(a40,i8)', ' # Cummulative number of vertices found:', nv_tot
    end if
  end if

  !-----------------------!
  !                       !
  !   Compress vertices   !
  !                       !
  !-----------------------!
  call Front % Compress_Vertices(verbose)

  !---------------------------------!
  !                                 !
  !   Find and check connectivity   !
  !                                 !
  !---------------------------------!
  call Front % Find_Sides(verbose)           ! Surf calls the same here
  call Front % Find_Front_Elements(verbose)  ! Surf calls Find_Surf_Elements
  call Front % Check_Elements(verbose)       ! Surf calls the same here

  !-----------------------------------------------!
  !   It used to find the nearest cell and node   !
  !   (But I deleted, I hope it was not needed)   !
  !-----------------------------------------------!
  n_verts_in_buffers = 0
  do v = 1, nv
    call Front % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers, locally=.true.)
    call Front % Vert(v) % Find_Nearest_Node()
  end do

  !--------------------------------------!
  !                                      !
  !   Calculate geometrical quantities   !
  !                                      !
  !--------------------------------------!
  call Front % Find_Vertex_Elements()
  call Front % Calculate_Element_Centroids()
  call Front % Calculate_Element_Normals(smooth)

  !-----------------------------------------------------------------!
  !                                                                 !
  !    Find cells at elements and intersect faces with elements     !
  !   (This is needed only for mass transfer problems with front)   !
  !                                                                 !
  !-----------------------------------------------------------------!
  call Front % Mark_Cells_And_Faces(sharp)

  call Cpu_Timer % Stop('Creating_Front_From_Vof_Function')

  return

  end subroutine
