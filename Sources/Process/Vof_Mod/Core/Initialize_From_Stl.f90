!==============================================================================!
  subroutine Initialize_From_Stl(Vof)
!------------------------------------------------------------------------------!
!   Initializes VOF function from an STL file                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NO    = 0
  integer, parameter :: YES   = 1
  logical, parameter :: DEBUG = .false.  ! if true, a lot of files are created
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer     :: Grid
  type(Stl_Type)               :: Stl
  type(Polyhedron_Type)        :: Pol(0:1)
  real,    allocatable         :: vof_body(:,:)
  real,    pointer, contiguous :: dis_nod_dom(:), dis_nod(:)
  real,    pointer, contiguous :: surf_dist_pos(:), surf_dist_neg(:)
  real,    pointer, contiguous :: surf_dist(:), node_dist(:)
  integer, pointer, contiguous :: dis_nod_int(:), cut_cel(:), glo(:)
  integer, pointer, contiguous :: set_cel(:), set_old(:)
  integer                      :: b, c, c1, c2, s, i_nod, j_nod, i, j, k1, k2, m
  integer                      :: i_iso, i_ver, i_fac, fac, p, reg
  integer                      :: n_cut_facets, cut_facets(1024)
  real                         :: vol_1, vol_2, vol_3, vol_4, vol_5
  real                         :: cell_vol, cel0_vol, cel1_vol, d
  real                         :: p1(3), p2(3), p3(3), qc(3), qn(3)
  real                         :: f(3), n(3)
  integer                      :: ij_cut(MAX_ISOAP_VERTS, MAX_ISOAP_VERTS)
  integer                      :: new_faces_n_nodes, cnt_p, cnt_m
  integer                      :: new_faces_n(MAX_ISOAP_VERTS)
  logical                      :: flooding
!==============================================================================!

  ! Take alias(es)
  Grid => Vof % pnt_grid
  glo  => Grid % Comm % cell_glo

  !-----------------------!
  !                       !
  !   Read the STL file   !
  !                       !
  !-----------------------!
  call Stl % Create_From_File(Vof % name_stl)

  allocate(vof_body(-Grid % n_bnd_cells:Grid % n_cells, Stl % n_boddies))

  !--------------------------------!
  !                                !
  !   Browse through STL boddies   !
  !                                !
  !--------------------------------!
  do b = 1, Stl % n_boddies

    ! In between boddies so that they reset all to zero
    call Work % Connect_Int_Cell(cut_cel, set_cel, set_old)
    call Work % Connect_Int_Node(dis_nod_int)
    call Work % Connect_Real_Node(dis_nod, dis_nod_dom)
    call Work % Connect_Real_Cell(surf_dist, surf_dist_pos, surf_dist_neg)
    call Work % Connect_Real_Cell(node_dist)

    if(Stl % n_boddies > 1 .and. First_Proc()) then
      print '(a,i3)', ' # Processing body ', b
    end if

    if(First_Proc()) then
      print '(a)', ' # Calculating distance from the STL interface'
    end if
    do c = Cells_In_Domain_And_Buffers()
      surf_dist_pos(c) = 0.0
      surf_dist_neg(c) = 0.0
      node_dist(c) = -HUGE
    end do

    do c = Cells_In_Domain()
      if(First_Proc()) then
        write(*,'(a2,f5.0,a14,a1)', advance='no') ' #',  &
             (100.*real(c)/real(Grid % n_cells)), ' % complete...', achar(13)
        flush(6)
      end if

      ! Node coordinates
      qc(1) = Grid % xc(c)
      qc(2) = Grid % yc(c)
      qc(3) = Grid % zc(c)

      ! Orthogonal distance to facets is computed in three steps:
      ! 1. find positive distance, this
      !    would be cells outside STL
      do fac = 1, Stl % N_Facets()
        if(Stl % Facet_Body(fac) .eq. b) then
          f = Stl % Facet_Coords(fac)
          n = Stl % Facet_Normal(fac)
          d = dot_product(qc-f, n)
          if(d > 0.0) then
            surf_dist_pos(c) = max(surf_dist_pos(c), d)
          end if
        end if  ! body
      end do
      ! 2. for cells which are left unmarked, hence
      !    inside the STL find negative distance
      if(surf_dist_pos(c) .eq. 0.0) then
        surf_dist_neg(c) = -HUGE
        do fac = 1, Stl % N_Facets()
          if(Stl % Facet_Body(fac) .eq. b) then
            f = Stl % Facet_Coords(fac)
            n = Stl % Facet_Normal(fac)
            d = dot_product(qc-f, n)
            if(d < 0.0) then
              surf_dist_neg(c) = max(surf_dist_neg(c), d)
            end if
          end if  ! body
        end do
      end if
      ! 3. work out the total distance
      surf_dist(c) = max(surf_dist_pos(c), abs(surf_dist_neg(c)))

      ! Cells internal dimensions
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        i = Grid % cells_n(i_nod, c)
        node_dist(c) = max(node_dist(c), Math % Distance_Squared(  &
          Grid % xc(c), Grid % yc(c), Grid % zc(c),                &
          Grid % xn(i), Grid % yn(i), Grid % zn(i)))
      end do
    end do
    call Grid % Exchange_Cells_Real(node_dist)
    call Grid % Exchange_Cells_Real(surf_dist)
    call Grid % Exchange_Cells_Real(surf_dist_pos)
    call Grid % Exchange_Cells_Real(surf_dist_neg)

    do c = Cells_In_Domain_And_Buffers()
      node_dist(c) = sqrt(node_dist(c))
    end do
    if(DEBUG) then
      call Grid % Save_Debug_Vtu(append="surf_dist_pos",       &
                                 scalar_cell=surf_dist_pos,    &
                                 scalar_name="surf_dist_pos")

      call Grid % Save_Debug_Vtu(append="surf_dist_neg",       &
                                 scalar_cell=surf_dist_neg,    &
                                 scalar_name="surf_dist_neg")

      call Grid % Save_Debug_Vtu(append="surf_dist",       &
                                 scalar_cell=surf_dist,    &
                                 scalar_name="surf_dist")

      call Grid % Save_Debug_Vtu(append="node_dist",          &
                                 scalar_cell=node_dist,       &
                                 scalar_name="node_dist")
    end if

    !----------------------------------------------------!
    !                                                    !
    !   Browse through cells to find out which are cut   !
    !   If a cell is cut, calculate distance from the    !
    !   STL file for all of its nodes.  Then normalize   !
    !   the node distance for all nodes in the grid      !
    !                                                    !
    !----------------------------------------------------!
    if(First_Proc()) print '(a)', ' # Searching for cells cut by the STL facets'
    do c = Cells_In_Domain_And_Buffers()

      ! Fetch cell coordinates
      qc(1) = Grid % xc(c);  qc(2) = Grid % yc(c);  qc(3) = Grid % zc(c)

      if(surf_dist(c) < 2.0 * node_dist(c)) then

        n_cut_facets  = 0  ! how many facets cut in this cell
        cut_facets(:) = 0
        cut_cel(c)    = NO

        if(First_Proc()) then
          write(*,'(a2,f5.0,a14,a1)', advance='no') ' #',  &
               (100. * real(c)/real(Grid % n_cells)), ' % complete...', achar(13)
          flush(6)
        end if

        !---------------------------------------------------------------!
        !   Browse through cells' nodes to find out if cells were cut   !
        !   (Not really all the cuts are detected in this way.  But     !
        !    it is OK, it misses only indents into edges)               !
        !---------------------------------------------------------------!
        do i_nod = 1, abs(Grid % cells_n_nodes(c))
          i = Grid % cells_n(i_nod, c)

          ! Fetch node coordinates
          qn(1) = Grid % xn(i);  qn(2) = Grid % yn(i);  qn(3) = Grid % zn(i)

          do fac = 1, Stl % N_Facets()
            if(Stl % Facet_Body(fac) .eq. b) then

              ! STL vertex coordinates (this is three times slower)
              p1(1:3) = Stl % Facets_Vert_Coords(fac, 1)
              p2(1:3) = Stl % Facets_Vert_Coords(fac, 2)
              p3(1:3) = Stl % Facets_Vert_Coords(fac, 3)

              ! Do i and j cross STL facet?
              vol_1 = dot_product(Math % Cross_Product(p1-qc, p2-qc), p3-qc)
              vol_2 = dot_product(Math % Cross_Product(p1-qn, p2-qn), p3-qn)
              ! vol_1 and vol_2 have different signs
              if(vol_1 * vol_2 < 0.0) then
                vol_3 = dot_product(Math % Cross_Product(qn-qc, p1-qc), p2-qc)
                vol_4 = dot_product(Math % Cross_Product(qn-qc, p2-qc), p3-qc)
                ! vol_3 and vol_3 have the same sign
                if(vol_3 * vol_4 > 0.0) then
                  vol_5 = dot_product(Math % Cross_Product(qn-qc, p3-qc), p1-qc)
                  ! vol_3, vol_4 and vol_5 have the same sign
                  if( (vol_3 < 0.0 .and. vol_4 < 0.0 .and. vol_5 < 0.0) .or.  &
                      (vol_3 > 0.0 .and. vol_4 > 0.0 .and. vol_5 > 0.0) ) then
                    n_cut_facets = n_cut_facets + 1
                    cut_facets(n_cut_facets) = fac
                    cut_cel(c) = YES
                  end if
                end if  ! vol_3 and vol_4 have the same sign
              end if  ! vol_1 and vol_2 have different signs
            end if  ! body

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
            end do  ! i_fac
          end do    ! i_nod

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

        end if  ! n_cut_facets > 0
      end if    ! surf_dist(c) < 2.0 * node_dist(c)
    end do      ! through cells

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
    if(DEBUG) then
      call Grid % Save_Debug_Vtu(append="dis_nod_dom",       &
                                 scalar_node=dis_nod_dom,    &
                                 scalar_name="dis_nod_dom")
      call Grid % Save_Debug_Vtu(append="cut_cell",             &
                                 scalar_cell=real(cut_cel,RP),  &
                                 scalar_name="cut_cel")
    end if

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
    if(DEBUG) then
      call Grid % Save_Debug_Vtu(append="dis_nod_int",              &
                                 scalar_node=real(dis_nod_int,RP),  &
                                 scalar_name="dis_nod_int")
    end if

    !---------------------------------------!
    !   The actual flood fill starts here   !
    !---------------------------------------!
    if(First_Proc()) write(*, '(a)', advance='no') ' # Flooding ...'
    m = 0
  1 continue
    m = m + 1
    flooding = .false.

    do c = Cells_In_Domain_And_Buffers()
      if( set_cel(c) .eq. NO ) then

        flooding = .true.     ! flooding is still going on
        vof_body(c, b) = 0.0  ! initialize to zero

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
          vof_body(c, b) = 1.0
          set_cel(c)     = YES  ! now it is set
          do i_nod = 1, abs(Grid % cells_n_nodes(c))
            i = Grid % cells_n(i_nod, c)
            dis_nod_int(i) = +1
          end do

        ! If there are more negative cell, set VOF to 0
        else if(cnt_p < cnt_m) then
          vof_body(c, b) = 0.0  ! probably not needed, it was set to 0 above
          set_cel(c)     = YES  ! now it is set
          do i_nod = 1, abs(Grid % cells_n_nodes(c))
            i = Grid % cells_n(i_nod, c)
            dis_nod_int(i) = -1
          end do

        ! Throw a warning if you come to this
        else   ! cnt_p == cnt_m
          if(cnt_p .gt. 0 .and. cnt_m .gt. 0) then
          ! PRINT *, __FILE__, __LINE__, DENTED CELL: ', C, CNT_P, CNT_M
          end if
        end if

      end if
    end do

    ! Store the old (from single processor) set_cel
    set_old(1:Grid % n_cells) = set_cel(1:Grid % n_cells)

    ! Take the set_cel and vof_bod(:,b) from other processors
    call Grid % Exchange_Cells_Int (set_cel)
    call Grid % Exchange_Cells_Real(vof_body(:,b))

    ! Check if the cell has been set in another processor
    do c = Cells_In_Domain_And_Buffers()
      if( set_cel(c) .eq. YES .and. set_old(c) .eq. NO ) then
        if(vof_body(c, b) < MICRO) then  ! the cell was filled with zero
          do i_nod = 1, abs(Grid % cells_n_nodes(c))
            i = Grid % cells_n(i_nod, c)
            dis_nod_int(i) = -1
          end do
        end if
        if(vof_body(c, b) > 1.0 - MICRO) then  ! the cell was filled with one
          do i_nod = 1, abs(Grid % cells_n_nodes(c))
            i = Grid % cells_n(i_nod, c)
            dis_nod_int(i) = +1
          end do
        end if
      end if
    end do

    ! Flood fill still going, go back
    call Global % Lor_Log(flooding)
    if(flooding) goto 1
    if(First_Proc()) print '(a)', ' done!'

    if(DEBUG) then
      call Grid % Save_Debug_Vtu(append="set_cel",                &
                                 scalar_cell=real(set_cel, RP),   &
                                 scalar_name="set_cel")
      call Grid % Save_Debug_Vtu(append="vof_tent",          &
                                 scalar_cell=vof_body(:,b),  &
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
    do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
      if( cut_cel(c) .eq. YES ) then

        !------------------------!
        !   Extract polyhedron   !
        !------------------------!
        call Polyhedron % Extract_From_Grid(Grid, c, dis_nod_dom)
        if(DEBUG) then
          call Polyhedron % Plot_Polyhedron_Vtk("dis-nod-dom", glo(c))
        end if

        !------------------------------!
        !   Call the Isoap algorithm   !
        !------------------------------!
        call Isoap % Extract_Iso_Polygons(Grid, c, dis_nod_dom)

        !---------------------------------------------------------------!
        !   It is indeed strange that here, during the initialization   !
        !    of VOF, some cells happen to have multiple iso-surfaces    !
        !---------------------------------------------------------------!
        if(Iso_Polygons % n_polys > 1) then
          print '(4(a,i8))',  __FILE__,        __LINE__,     &
                             ' # check cell ', c,            &
                             ' in processor ', This_Proc(),  &
                             ' global cell ',  glo(c)
          call Polyhedron % Plot_Polyhedron_Vtk("check-cell", glo(c))
          call Iso_Polygons % Plot_Iso_polygons_Vtk("check-iso", glo(c))
        end if

        !--------------------------------!
        !   Add new vertices for faces   !
        !--------------------------------!
        ij_cut(:,:) = 0
        do s = 1, Polyhedron % n_faces

          new_faces_n_nodes = 0  ! initalize number of nodes in new face
          new_faces_n(:)    = 0  ! initialize new face's node list to zero

          do i_nod = 1, Polyhedron % faces_n_nodes(s)
            j_nod = i_nod + 1
            if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1

            ! Two nodes at polyhedron
            i = Polyhedron % faces_n(s, i_nod)
            j = Polyhedron % faces_n(s, j_nod)

            ! Plain copy of the old node (i) to new face's nodes
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n(new_faces_n_nodes) = i  ! just copy "i"

            do i_iso = 1, Iso_Polygons % n_polys
              do i_ver = 1, Iso_Polygons % polys_n_verts(i_iso)
                m = Iso_Polygons % polys_v(i_iso, i_ver)

                k1 = Iso_Polygons % b_node_1(m)
                k2 = Iso_Polygons % b_node_2(m)

                if(k1 == i .and. k2 == j .or.  &
                   k1 == j .and. k2 == i) then

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
            call Polyhedron % Plot_Polyhedron_Vtk("cell", glo(c))
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
              if( p .eq. 0 .and. Pol(p) % phi(i) < 0.5 + MICRO .or.  &
                  p .eq. 1 .and. Pol(p) % phi(i) > 0.5 - MICRO) then
                new_faces_n_nodes = new_faces_n_nodes + 1
                new_faces_n  (new_faces_n_nodes) = i  ! just copy "i"
              end if
            end do  ! through nodes of the face

            ! Now when all the nodes have been browsed, reform the
            ! face, I mean overwrite the old one with the new one
            Pol(p) % faces_n_nodes(s) = new_faces_n_nodes ! copy number of nodes
            Pol(p) % faces_n(s,1:new_faces_n_nodes)  &    ! copy list of nodes
                 = new_faces_n(1:new_faces_n_nodes)
          end do    ! through polyhedron faces
          ! What if there is a face without any nodes between 0 and 0.5+MICRO?
          ! Well, faces_n_nodes will be zero for that face

          if(DEBUG) then
            if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0-hollow",glo(c))
            if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1-hollow",glo(c))
          end if

          ! Try to add the missing face
          new_faces_n_nodes = 0
          new_faces_n(:)    = 0
          do i_iso = 1, Iso_Polygons % n_polys
            do i_ver = 1, Iso_Polygons % polys_n_verts(i_iso)
              m = Iso_Polygons % polys_v(i_iso, i_ver)
              do i = 1, Pol(p) % n_nodes
                if(Math % Approx_Real(Iso_Polygons % verts_xyz(m,1),        &
                                            Pol(p) % nodes_xyz(i,1)) .and.  &
                   Math % Approx_Real(Iso_Polygons % verts_xyz(m,2),        &
                                            Pol(p) % nodes_xyz(i,2)) .and.  &
                   Math % Approx_Real(Iso_Polygons % verts_xyz(m,3),        &
                                            Pol(p) % nodes_xyz(i,3)) ) then
                  new_faces_n_nodes = new_faces_n_nodes + 1
                  new_faces_n(new_faces_n_nodes) = i
                end if  ! node/vert match
              end do  ! i
            end do    ! i_ver
          end do      ! i_iso
          Pol(p) % n_faces = Pol(p) % n_faces + 1        ! one more face
          s = Pol(p) % n_faces                           ! to shorten the syntax
          Pol(p) % faces_n_nodes(s) = new_faces_n_nodes  ! copy number of nodes
          Pol(p) % faces_n(s,1:new_faces_n_nodes)  &     ! copy list of nodes
             = new_faces_n(1:new_faces_n_nodes)

          if(DEBUG) then
            if(p .eq. 0) call Pol(0) % Plot_Polyhedron_Vtk("pol0-full", glo(c))
            if(p .eq. 1) call Pol(1) % Plot_Polyhedron_Vtk("pol1-full", glo(c))
          end if
        end do  ! through p

        ! Calculate volume of the Pol(0)
        call Pol(0)     % Calculate_Cell_Volume(cel0_vol)
        call Pol(1)     % Calculate_Cell_Volume(cel1_vol)
        call Polyhedron % Calculate_Cell_Volume(cell_vol)
        if(cel0_vol > cel0_vol) then
          vof_body(c, b) = (cell_vol - cel0_vol) / (cell_vol)
        else
          vof_body(c, b) = cel1_vol / (cell_vol)
        end if


      end if  ! cut_cel(c) .eq. YES
    end do

    if(b .eq. 1) then
      do c = 1, Grid % n_cells
        Vof % fun % n(c) = vof_body(c, b)
      end do
    else
      do c = 1, Grid % n_cells
        Vof % fun % n(c) = Vof % fun % n(c) * vof_body(c, b)
      end do
    end if

    ! In between boddies so that they reset all to zero
    call Work % Disconnect_Int_Cell(cut_cel, set_cel, set_old)
    call Work % Disconnect_Int_Node(dis_nod_int)
    call Work % Disconnect_Real_Node(dis_nod, dis_nod_dom)
    call Work % Disconnect_Real_Cell(surf_dist, surf_dist_pos, surf_dist_neg)
    call Work % Disconnect_Real_Cell(node_dist)

  end do  ! through boddies
  call Grid % Exchange_Cells_Real(Vof % fun % n)

  !-----------------------------!
  !   Set boundary values too   !
  !-----------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. INFLOW) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Vof % fun % n(c2) = Vof % fun % n(c1)
      end do  ! faces
    end if    ! inflow
  end do      ! region

  !------------------------------------!
  !   Set the old time steps as well   !
  !------------------------------------!
  do reg = Boundary_Inside_And_Buffer_Regions()
    do c = Cells_In_Region(reg)
      Vof % fun % o (c) = Vof % fun % n(c)
      Vof % fun % oo(c) = Vof % fun % n(c)
    end do
  end do

  if(DEBUG) then
    call Grid % Save_Debug_Vtu(append="vof_final",         &
                               scalar_cell=Vof % fun % n,  &
                               scalar_name="vof_final")
  end if

  !--------------------!
  !                    !
  !   Smooth the VOF   !
  !                    !
  !--------------------!
  call Vof % Smooth_Vof_And_Compute_Surface_Normals()

  !----------------------!
  !                      !
  !   Initialize front   !
  !                      !
  !----------------------!
  if(Vof % track_front) then
    call Vof % Front % Place_Front_At_Value(Vof % fun, .true.)
    call Vof % Front % Print_Front_Statistics()
  end if

  !---------------------------------!
  !                                 !
  !   It was initialized from STL   !
  !                                 !
  !---------------------------------!
  Vof % init_stl = .true.

  end subroutine
