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
  integer, parameter :: NO    = 0
  integer, parameter :: YES   = 1
  logical, parameter :: DEBUG = .true.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Stl_Type)               :: Stl
  type(Polyhedron_Type)        :: Pol(0:1)
  real,    pointer, contiguous :: dis_nod_dom(:), dis_nod(:)
  integer, pointer, contiguous :: dis_nod_int(:), cut_cel(:), set_cel(:)
  integer                      :: c, s, i_nod, j_nod, i, j, k1, k2, m
  integer                      :: i_iso, i_ver, i_fac, fac, p
  integer                      :: n_cut_facets, cut_facets(1024)
  real                         :: vol_1, vol_2, vol_3, vol_4, vol_5
  real                         :: cell_vol, cel0_vol, cel1_vol
  real                         :: p1(3), p2(3), p3(3), qi(3), qj(3), qn(3)
  real                         :: f(3), n(3)
  integer                      :: ij_cut(MAX_ISOAP_VERTS, MAX_ISOAP_VERTS)
  integer                      :: new_faces_n_nodes, cnt_p, cnt_m
  integer                      :: new_faces_n(MAX_ISOAP_VERTS)
  logical                      :: flooding
!==============================================================================!

  call Work % Connect_Int_Cell(cut_cel)
  call Work % Connect_Int_Cell(set_cel)
  call Work % Connect_Int_Node(dis_nod_int)
  call Work % Connect_Real_Node(dis_nod)
  call Work % Connect_Real_Node(dis_nod_dom)

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !-----------------------!
  !                       !
  !   Read the STL file   !
  !                       !
  !-----------------------!
  call Stl % Create_From_File("2_spheres.stl")

  !----------------------------------------------------!
  !                                                    !
  !   Browse through cells to find out which are cut   !
  !   If a cell is cut, calculate distance from the    !
  !   STL file for all of its nodes.  Then normalize   !
  !   the node distance for all nodes in the grid      !
  !                                                    !
  !----------------------------------------------------!
  do c = 1, Grid % n_cells

    n_cut_facets  = 0  ! how many facets cut in this cell
    cut_facets(:) = 0
    cut_cel(c)    = NO

    write(*,'(a2,f5.0,a14,a1)', advance='no')  &
      ' #', (100. * real(c)/real(Grid % n_cells)), ' % complete...', achar(13)

    !---------------------------------------------------------------!
    !   Browse through cells' nodes to find out if cells were cut   !
    !   (Not really all the cuts are detected in this way.  But     !
    !    it is OK, it misses only indents into edges)               !
    !---------------------------------------------------------------!
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      i = Grid % cells_n(i_nod, c)

      qi(1) = Grid % xc(c)
      qi(2) = Grid % yc(c)
      qi(3) = Grid % zc(c)
      qj(1) = Grid % xn(i)
      qj(2) = Grid % yn(i)
      qj(3) = Grid % zn(i)

      do fac = 1, Stl % N_Facets()

        ! STL vertex coordinates (this is three times slower)
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
        qn(1) = Grid % xn(i)
        qn(2) = Grid % yn(i)
        qn(3) = Grid % zn(i)

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
        dis_nod_int(i) = dis_nod_int(i) + 1
      end do

    end if
  end do    ! through cells

  !--------------------------------------------------!
  !   Equalize distance at nodes for entire domain   !
  !   (Note that this will spread distances a bit    !
  !    even to the cells which are not cut.)         !
  !--------------------------------------------------!
  do i = 1, Grid % n_nodes
    if(dis_nod_int(i) > 0) then
      dis_nod_dom(i) = dis_nod_dom(i) / real(dis_nod_int(i))
    else
      dis_nod_dom(i) = 0.0
    end if
  end do

  ! See what you have at this moment
  call Grid % Save_Debug_Vtu(append="dis_nod_dom",       &
                             scalar_node=dis_nod_dom,    &
                             scalar_name="dis_nod_dom")
  call Grid % Save_Debug_Vtu(append="cut_cell",             &
                             scalar_cell=real(cut_cel,RP),  &
                             scalar_name="cut_cel")

  !-------------------------------------!
  !                                     !
  !   Find vof in cells by flood fill   !
  !                                     !
  !-------------------------------------!

  ! Set set_cel to be the same as cut_cel as ...
  ! ... cells which are not cut are also not set
  set_cel(1:Grid % n_cells) = cut_cel(1:Grid % n_cells)

  ! Turn distances to integers
  do i = 1, Grid % n_nodes
    if(Math % Approx_Real(dis_nod_dom(i), 0.0)) then
      dis_nod_int(i) = 0
    else
      if(dis_nod_dom(i) > 0.0) then
        dis_nod_int(i) = 1
      else
        dis_nod_int(i) = -1
      end if
    end if
  end do
  call Grid % Save_Debug_Vtu(append="dis_nod_int",              &
                             scalar_node=real(dis_nod_int,RP),  &
                             scalar_name="dis_nod_int")

  !---------------------------------------!
  !   The actual flood fill starts here   !
  !---------------------------------------!
  print '(a)', ' # Flooding ...'
  m = 0
1 continue
  m = m + 1
  flooding = .false.
    write(*,'(a2,f5.0,a14,a1)', advance='no')  &
      ' #', (100. * real(sum(set_cel))/real(Grid % n_cells)), ' % complete...', achar(13)

  do c = 1, Grid % n_cells
    if( set_cel(c) .eq. NO ) then

      flooding = .true.      ! flooding is still going on
      Vof % fun % n(c) = 0.0  ! initialize to zero

      ! Count nodes with integer distance -1 and +1
      cnt_m = 0
      cnt_p = 0
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        i = Grid % cells_n(i_nod, c)
        if(dis_nod_int(i) .eq. +1) cnt_p = cnt_p + 1
        if(dis_nod_int(i) .eq. -1) cnt_m = cnt_m + 1
      end do

      ! If there are more positive cell, set VOF to 1
      if(cnt_p > cnt_m) then
        Vof % fun % n(c) = 1.0
        set_cel(c)      = YES  ! now it is set
        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          i = Grid % cells_n(i_nod, c)
          dis_nod_int(i) = +1
        end do

      ! If there are more negative cell, set VOF to 0
      else if(cnt_p < cnt_m) then
        Vof % fun % n(c) = 0.0  ! probably not needed, it was set to 0 above
        set_cel(c)       = YES  ! now it is set
        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          i = Grid % cells_n(i_nod, c)
          dis_nod_int(i) = -1
        end do

      ! Throw an error if you come to this
      else   ! cnt_p == cnt_m
        if(cnt_p .gt. 0 .and. cnt_m .gt. 0) then
        PRINT *, 'HOW ON EARTH?  CELL: ', C, CNT_P, CNT_M
        end if
      end if

    end if
  end do

  ! Flood fill still going, go back
  if(flooding) goto 1

  ! Exchange is needed here
  if(DEBUG) then
    call Grid % Save_Debug_Vtu(append="set_cel",                &
                               scalar_cell=real(set_cel, RP),   &
                               scalar_name="set_cel")
    call Grid % Save_Debug_Vtu(append="vof_tent",          &
                               scalar_cell=Vof % fun % n,  &
                               scalar_name="vof_tent")
  end if

  !------------------------------------!
  !                                    !
  !   Calculate VOF in all cut cells   !
  !                                    !
  !------------------------------------!

  ! Set dis_nod_dom to usable values; values needed by Isoap
  do i = 1, Grid % n_nodes
    dis_nod_dom(i) = dis_nod_dom(i) + 0.5
  end do

  !-----------------------------------!
  !   Browse through cut cells only   !
  !-----------------------------------!
  do c = 1, Grid % n_cells
    if( cut_cel(c) .eq. YES ) then

      !------------------------!
      !   Extract polyhedron   !
      !------------------------!
      call Polyhedron % Extract_From_Grid(Grid, c, dis_nod_dom)
      if(DEBUG) then
        call Polyhedron % Plot_Polyhedron_Vtk("dis_nod_dom", c)
      end if

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
        if(DEBUG) then
          call Polyhedron % Plot_Polyhedron_Vtk("geo", c)
        end if

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

        if(DEBUG) then
          if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0_hollow", c)
          if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1_hollow", c)
        end if

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

        if(DEBUG) then
          if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0_full", c)
          if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1_full", c)
        end if
      end do  ! through p

      ! Calculate volume of the Pol(0)
      call Pol(0)     % Calculate_Cell_Volume(cel0_vol)
      call Pol(1)     % Calculate_Cell_Volume(cel1_vol)
      call Polyhedron % Calculate_Cell_Volume(cell_vol)
      if(cel0_vol > cel0_vol) then
        Vof % fun % n(c) = (cell_vol - cel0_vol) / (cell_vol)
      else
        Vof % fun % n(c) = cel1_vol / (cell_vol)
      end if

    end if  ! cut_cel(c) .eq. YES
  end do

  ! Exchange is needed here
  call Grid % Save_Debug_Vtu(append="vof_final",         &
                             scalar_cell=Vof % fun % n,  &
                             scalar_name="vof_final")

  call Work % Disconnect_Int_Cell(cut_cel)
  call Work % Disconnect_Int_Cell(set_cel)
  call Work % Disconnect_Int_Node(dis_nod_int)
  call Work % Disconnect_Real_Node(dis_nod)
  call Work % Disconnect_Real_Node(dis_nod_dom)

  call Comm_Mod_End
  stop

  end subroutine
