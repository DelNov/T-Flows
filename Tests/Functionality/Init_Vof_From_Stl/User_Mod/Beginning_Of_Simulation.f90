!==============================================================================!
  real function Sgn_Volume(p1, p2, p3, p4)
!------------------------------------------------------------------------------!
!   Returns the signed volume of tethraedra spanned with nodes 1 to 4          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real, intent(in) :: p1(3), p2(3), p3(3), p4(3)
!==============================================================================!
  Sgn_Volume = (   (  (p2(2)-p1(2))*(p3(3)-p1(3))                      &
                    - (p3(2)-p1(2))*(p2(3)-p1(3)) ) * (p4(1)-p1(1)) +  &
                   (  (p3(1)-p1(1))*(p2(3)-p1(3))                      &
                    - (p2(1)-p1(1))*(p3(3)-p1(3)) ) * (p4(2)-p1(2)) +  &
                   (  (p2(1)-p1(1))*(p3(2)-p1(2))                      &
                    - (p3(1)-p1(1))*(p2(2)-p1(2)) ) * (p4(3)-p1(3)) ) / 6.0
  end function

!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NO  = 0
  integer, parameter :: YES = 1
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Stl_Type)               :: Stl
  type(Polyhedron_Type)        :: Pol(0:1)
  real,    pointer, contiguous :: vof_at_nodes(:), vof_in_cells(:)
  real,    pointer, contiguous :: dis_nod_dom(:), dis_nod(:)
  integer, pointer, contiguous :: dis_nod_cnt(:), cut_cel(:)
  integer                      :: c, s, i_nod, j_nod, i, j, k1, k2, m
  integer                      :: i_iso, i_ver, i_fac, fac, p, n_cut_facets
  integer                      :: cut_facets(1024)
  real                         :: vol_1, vol_2, vol_3, vol_4, vol_5
  real                         :: cell_vol, cel0_vol, cel1_vol
  real                         :: p1(3), p2(3), p3(3), qi(3), qj(3), qn(3)
  real                         :: f(3), n(3), l(3)
  integer                      :: ij_cut(MAX_ISOAP_VERTS, MAX_ISOAP_VERTS)
  integer                      :: new_faces_n_nodes
  integer                      :: new_faces_n(MAX_ISOAP_VERTS)
!==============================================================================!

  call Work % Connect_Int_Node(dis_nod_cnt)
  call Work % Connect_Real_Node(vof_at_nodes, dis_nod, dis_nod_dom)
  call Work % Connect_Int_Cell(cut_cel)
  call Work % Connect_Real_Cell(vof_in_cells)

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !-----------------------!
  !   Read the STL file   !
  !-----------------------!
  call Stl % Create_From_File("2_spheres.stl")

  !-----------------------!
  !   Find vof in nodes   !
  !-----------------------!
  do i = 1, Grid % n_nodes

    vof_at_nodes(i) = 0  ! assume node is in the stl

    do fac = 1, Stl % N_Facets()

      f = Stl % Facet_Coords(fac)
      n = Stl % Facet_Normal(fac)

      ! Vector connecting facet centroid with the node
      qn(1:3) = (/Grid % xn(i)-f(1), Grid % yn(i)-f(2), Grid % zn(i)-f(3)/)

      ! First time this product is positive, node is outside of STL
      if(dot_product(qn, n) > 0) then
        vof_at_nodes(i) = vof_at_nodes(i) + 1.0
      end if
    end do

    vof_at_nodes(i) = vof_at_nodes(i) / real(Stl % N_Facets())
    if(vof_at_nodes(i) < 0.01) then
      vof_at_nodes(i) = 0.0
    else
      vof_at_nodes(i) = 1.0
    end if
  end do
  call Grid % Save_Debug_Vtu(append="vof_at_nodes",       &
                             scalar_node=vof_at_nodes,    &
                             scalar_name="vof_at_nodes")

  !-----------------------!
  !   Find vof in cells   !
  !-----------------------!
  do c = 1, Grid % n_cells
    vof_in_cells(c) = 0.0
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      i = Grid % cells_n(i_nod, c)
      vof_in_cells(c) = vof_in_cells(c) + vof_at_nodes(i)
    end do
    vof_in_cells(c) = vof_in_cells(c) / real(abs(Grid % cells_n_nodes(c)))
  end do
  ! Exchange is needed here
  call Grid % Save_Debug_Vtu(append="vof_in_cells_tentative",  &
                             scalar_cell=vof_in_cells,         &
                             scalar_name="vof_in_cells")

  !--------------------------!
  !   Browse through cells   !
  !--------------------------!
  do c = 1, Grid % n_cells

    n_cut_facets  = 0  ! how many facets cut in this cell
    cut_facets(:) = 0
    cut_cel(c) = NO

    !---------------------------------------------------------------!
    !   Browse through cells' nodes to find out if cells were cut   !
    !   (Not really all the cuts are detected in this way.  But     !
    !    it is OK, it misses only indents into edges)               !
    !---------------------------------------------------------------!
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      i = Grid % cells_n(i_nod, c)

      qi(1:3) = (/Grid % xc(c), Grid % yc(c), Grid % zc(c)/)
      qj(1:3) = (/Grid % xn(i), Grid % yn(i), Grid % zn(i)/)

      do fac = 1, Stl % N_Facets()

        ! STL vertex coordinates
        p1(1:3) = Stl % Facets_Vert_Coords(fac, 1)
        p2(1:3) = Stl % Facets_Vert_Coords(fac, 2)
        p3(1:3) = Stl % Facets_Vert_Coords(fac, 3)

        ! Do i and j cross STL facet?
        vol_1 = Sgn_Volume(qi, p1, p2, p3)
        vol_2 = Sgn_Volume(qj, p1, p2, p3)
        ! vol_1 and vol_2 have different signs
        if(vol_1 * vol_2 < 0.0) then
          vol_3 = Sgn_Volume(qi, qj, p1, p2)
          vol_4 = Sgn_Volume(qi, qj, p2, p3)
          ! vol_3 and vol_3 have the same sign
          if(vol_3 * vol_4 > 0.0) then
            vol_5 = Sgn_Volume(qi, qj, p3, p1)
            ! vol_3, vol_4 and vol_5 have the same sign
            if( (vol_3 < 0.0 .and. vol_4 < 0.0 .and. vol_5 < 0.0) .or.  &
                (vol_3 > 0.0 .and. vol_4 > 0.0 .and. vol_5 > 0.0) ) then
              n_cut_facets = n_cut_facets + 1
              cut_facets(n_cut_facets) = fac
              cut_cel(c) = YES
            end if
          end if  ! vol_3 and vol_4 have the same sign
        end if  ! vol_1 and vol_2 have different signs

      end do  ! STL facets
    end do    ! cells' nodes

    !----------------------------------------------!
    !   If cell was cut, calculate distance from   !
    !      each of its nodes to the interface      !
    !----------------------------------------------!
    if(n_cut_facets > 0) then  ! and cut_cel(c) .eq. YES

      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        i = Grid % cells_n(i_nod, c)

        dis_nod(i) = 0.0

        ! Node coordinates
        qn(1:3) = (/Grid % xn(i), Grid % yn(i), Grid % zn(i)/)

        ! Accumulate distance to the STL facets
        do i_fac = 1, n_cut_facets
          fac = cut_facets(i_fac)
          f   = Stl % Facet_Coords(fac)
          n   = Stl % Facet_Normal(fac)
          dis_nod(i) = dis_nod(i) + dot_product(qn-f, n)
        end do
      end do

      ! Average the distances from all STL facets
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        i = Grid % cells_n(i_nod, c)
        dis_nod(i) = dis_nod(i) / real(n_cut_facets)
      end do

      ! Add this local, cell-wise distances to that of the domain
      ! (This could be done together with the loop above this)
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        i = Grid % cells_n(i_nod, c)
        dis_nod_dom(i) = dis_nod_dom(i) + dis_nod(i)
        dis_nod_cnt(i) = dis_nod_cnt(i) + 1
      end do

    end if
  end do    ! through cells

  !------------------------------------------------------!
  !   Equalize distance at nodes for the entire domain   !
  !------------------------------------------------------!
  do i = 1, Grid % n_nodes
    if(dis_nod_cnt(i) > 0) then
      dis_nod_dom(i) = dis_nod_dom(i) / real(dis_nod_cnt(i))
    else
      dis_nod_dom(i) = 0.0
    end if
  end do

  ! See what you have at this moment
  call Grid % Save_Debug_Vtu(append="dis_nod_dom",       &
                             scalar_node=dis_nod_dom,    &
                             scalar_name="dis_nod_dom")

  ! Set dis_nod_dom to usable values
  do i = 1, Grid % n_nodes
    dis_nod_dom(i) = dis_nod_dom(i) + 0.5
  end do

  !------------------------------------!
  !                                    !
  !   Browse through cut cells again   !
  !                                    !
  !------------------------------------!
  do c = 1, Grid % n_cells
    if( cut_cel(c) .eq. YES ) then

      !------------------------!
      !   Extract polyhedron   !
      !------------------------!
      call Polyhedron % Extract_From_Grid(Grid, c, dis_nod_dom)
      call Polyhedron % Plot_Polyhedron_Vtk("dis_nod_dom", c)

      !------------------------------!
      !   Call the Isoap algorithm   !
      !------------------------------!
      call Isoap % Extract_Iso_Polygons(Grid, c, dis_nod_dom)
      call Iso_Polygons % Plot_Iso_Polygons_Vtk(c)

      !--------------------------------!
      !   Add new vertices for faces   !
      !--------------------------------!
      ij_cut(:,:) = 0
      do s = 1, Polyhedron % n_faces

        new_faces_n_nodes = 0  ! initalize number of nodes in new face
        new_faces_n(:)    = 0  ! initialize new face's node list to zero

        do i_nod = 1, Polyhedron % faces_n_nodes(s)
          j_nod = i_nod + 1; if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1

          ! Two nodes at polyhedron
          i = Polyhedron % faces_n(s, i_nod)
          j = Polyhedron % faces_n(s, j_nod)

          ! Plain copy of the old node (i) to new face's nodes
          new_faces_n_nodes = new_faces_n_nodes + 1
          new_faces_n(new_faces_n_nodes) = i  ! just copy "i"

          IF(ISO_POLYGONS % N_POLYS > 1) PRINT *, 'CHECK CELL: ', C
          do i_iso = 1, Iso_Polygons % n_polys
            do i_ver = 1, Iso_Polygons % polys_n_verts(i_iso)
              m = Iso_Polygons % polys_v(i_iso, i_ver)

              k1 = Iso_Polygons % b_node_1(m)
              k2 = Iso_Polygons % b_node_2(m)

              if(k1 == i .and. k2 == j .or.  &
                 k1 == j .and. k2 == i) then
                ! PRINT *, Iso_Polygons % verts_xyz(m, 1:3)

                ! Node hasn't been added yet, add it to the polygon
                if(ij_cut(j,i) .eq. 0) then
                  Polyhedron % n_nodes = Polyhedron % n_nodes + 1
                  Polyhedron % nodes_xyz(Polyhedron % n_nodes, 1:3)   &
                                 = Iso_Polygons % verts_xyz(m, 1:3)
                  Polyhedron % phi(Polyhedron % n_nodes) = 0.5

                  new_faces_n_nodes = new_faces_n_nodes + 1
                  new_faces_n(new_faces_n_nodes) = Polyhedron % n_nodes

                  ij_cut(i,j) = Polyhedron % n_nodes

                ! Node has already been added, don't add the node,
                ! but add it to the list of nodes in this face
                else
                  new_faces_n_nodes = new_faces_n_nodes + 1
                  new_faces_n(new_faces_n_nodes) = ij_cut(j,i)
                end if

              end if
            end do
          end do

        end do  ! i_nod

        ! Now when all the nodes have been browsed, reform the
        ! face, I mean overwrite the old one with the new one
        Polyhedron % faces_n_nodes(s) = new_faces_n_nodes
        Polyhedron % faces_n(s, 1:new_faces_n_nodes)  &
                  = new_faces_n(1:new_faces_n_nodes)
        call Polyhedron % Plot_Polyhedron_Vtk("geo", c)

      end do  ! s

      !------------------------------!
      !   Create Pol(0) and Pol(1)   !
      !------------------------------!
      do p = 0, 1
        call Pol(p) % Create_From_Polyhedron(Polyhedron)

        do s = 1, Pol(p) % n_faces
          new_faces_n_nodes = 0
          new_faces_n(:)    = 0

          ! Browse nodes in circular direction
          do i_nod = 1, Pol(p) % faces_n_nodes(s)
            i = Pol(p) % faces_n(s, i_nod)
            if( p .eq. 0 .and. Pol(p) % phi(i) < 0.5+MICRO .or.  &
                p .eq. 1 .and. Pol(p) % phi(i) > 0.5-MICRO) then
              new_faces_n_nodes = new_faces_n_nodes + 1
              new_faces_n  (new_faces_n_nodes) = i  ! just copy "i"
            end if
          end do  ! through nodes of the face

          ! Now when all the nodes have been browsed, reform the
          ! face, I mean overwrite the old one with the new one
          Pol(p) % faces_n_nodes(s) = new_faces_n_nodes ! copy number of nodes
          Pol(p) % faces_n(s,1:new_faces_n_nodes)  &    ! copy the list of nodes
               = new_faces_n(1:new_faces_n_nodes)
        end do    ! through polyhedron faces
        ! What if there is a face without any nodes between 0 and 0.5+MICRO?
        ! Well, faces_n_nodes will be zero for that face

        if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0_hollow", c)
        if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1_hollow", c)

        ! Try to add the missing face
        new_faces_n_nodes = 0
        new_faces_n(:)    = 0
        do i_iso = 1, Iso_Polygons % n_polys
          do i_ver = 1, Iso_Polygons % polys_n_verts(i_iso)
            m = Iso_Polygons % polys_v(i_iso, i_ver)
            do i = 1, Pol(p) % n_nodes
              if(Math % Approx_Three_Reals(                   &
                          Iso_Polygons % verts_xyz(m, 1:3),   &
                            Pol(p) % nodes_xyz(i, 1:3))) then
                new_faces_n_nodes = new_faces_n_nodes + 1
                new_faces_n(new_faces_n_nodes) = i
              end if  ! node/vert match
            end do  ! i
          end do    ! i_ver
        end do      ! i_iso
        Pol(p) % n_faces = Pol(p) % n_faces + 1        ! one more face
        s = Pol(p) % n_faces                           ! to shorten the syntax
        Pol(p) % faces_n_nodes(s) = new_faces_n_nodes  ! copy number of nodes
        Pol(p) % faces_n(s,1:new_faces_n_nodes)  &     ! copy the list of nodes
           = new_faces_n(1:new_faces_n_nodes)

        if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0_full", c)
        if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1_full", c)
      end do  ! through p

      ! Calculate volume of the Pol(0)
      call Pol(0)     % Calculate_Cell_Volume(cel0_vol)
      call Pol(1)     % Calculate_Cell_Volume(cel1_vol)
      call Polyhedron % Calculate_Cell_Volume(cell_vol)
      if(cel0_vol > cel0_vol) then
        vof_in_cells(c) = (cell_vol - cel0_vol) / (cell_vol)
      else
        vof_in_cells(c) = cel1_vol / (cell_vol)
      end if

    end if  ! cut_cel(c) .eq. YES
  end do

  ! Exchange is needed here
  call Grid % Save_Debug_Vtu(append="vof_in_cells_final",  &
                             scalar_cell=vof_in_cells,     &
                             scalar_name="vof_in_cells")

  call Work % Disconnect_Int_Node(dis_nod_cnt)
  call Work % Disconnect_Real_Node(vof_at_nodes, dis_nod, dis_nod_dom)
  call Work % Disconnect_Real_Cell(vof_in_cells)

  call Comm_Mod_End
  stop

  end subroutine
