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
  integer, parameter :: UNDEFINED = 15
  integer, parameter :: LOW       =  0
  integer, parameter :: HIGH      = 10
  integer, parameter :: MEDIUM    =  5
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Stl_Type)           :: Sphere
  type(Polyhedron_Type)    :: Pol0
  type(Polyhedron_Type)    :: Pol1
  real, pointer, contiguous :: vof_nod(:), vof_cel(:)
  integer                  :: cut_count, new_faces_n_nodes, n_combo, cc, next_n
  real                     :: tot_vol, cel_vol, di
  real                     :: vol_1, vol_2, vol_3, vol_4, vol_5
  real                     :: p1(3), p2(3), p3(3), qi(3), qj(3), qn(3)
  real                     :: f(3), n(3), l(3), xyz_new(3)
  integer                  :: c, s, fac, i, j, i_nod, j_nod, i_fac, run
  logical                  :: ij_cut_flag, has_one_point_five
  integer                  :: ij_cut(MAX_ISOAP_VERTS, MAX_ISOAP_VERTS)
  integer                  :: new_faces_n(MAX_ISOAP_VERTS)
  integer                  :: combos(MAX_ISOAP_FACES, 2)
!==============================================================================!

  call Work % Connect_Real_Node(vof_nod)
  call Work % Connect_Real_Cell(vof_cel)

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !-------------------------------------------------------------!
  !   Calculate and report integral volume in traditional way   !
  !-------------------------------------------------------------!
  tot_vol = 0.0
  do c = 1, Grid % n_cells
    tot_vol = tot_vol + Grid % vol(c)
  end do
  print '(a, f12.5)', ' # Total volume from Grid is: ', tot_vol

  do c = 1, Grid % n_cells
    call Polyhedron % Extract_From_Grid(Grid, c)
    ! if(mod(c,100) .eq. 0) call Polyhedron % Plot_Polyhedron_Vtk(c)
  end do

  !---------------------------------------!
  !   Calculate individual cell volumes   !
  !---------------------------------------!
  tot_vol = 0.0
  do c = 1, Grid % n_cells
    call Polyhedron % Extract_From_Grid(Grid, c)
    call Polyhedron % Calculate_Cell_Volume(cel_vol)
    tot_vol = tot_vol + cel_vol
  end do
  print '(a, f12.5)', ' # Total volume from Polyhedron is: ', tot_vol

  !------------------------------!
  !                              !
  !   Now comes the real stuff   !
  !                              !
  !------------------------------!
  call Sphere % Create_From_File("sphere.stl")

  !-----------------------!
  !   Find vof in nodes   !
  !-----------------------!
  do i = 1, Grid % n_nodes

    vof_nod(i) = 0  ! assume node is in the stl

    do fac = 1, Sphere % n_facets

      ! STL vertex coordinates
      p1(1)=Sphere % x(1,fac); p1(2)=Sphere % y(1,fac); p1(3)=Sphere % z(1,fac)
      p2(1)=Sphere % x(2,fac); p2(2)=Sphere % y(2,fac); p2(3)=Sphere % z(2,fac)
      p3(1)=Sphere % x(3,fac); p3(2)=Sphere % y(3,fac); p3(3)=Sphere % z(3,fac)

      ! STL facet centroid
      f = (p1+p2+p3)/3.0

      ! STL facet normal
      n(1)=Sphere % nx(fac);  n(2)=Sphere % ny(fac);  n(3)=Sphere % nz(fac)

      if(dot_product(f, n) < 0.0) then
        PRINT *, 'OUCH, THE ORIENTATION oF THE NORMAL IS WRONG'
      end if

      ! Vector connecting facet centroid with the node
      qn(1) = Grid % xn(i) - f(1)
      qn(2) = Grid % yn(i) - f(2)
      qn(3) = Grid % zn(i) - f(3)

      ! First time this product is positive, node is outside of STL
      if(dot_product(qn, n) > 0) then
        vof_nod(i) = vof_nod(i) + 1.0
      end if
    end do

    vof_nod(i) = vof_nod(i) / real(Sphere % n_facets)
    if(vof_nod(i) < 0.01) then
      vof_nod(i) = 0.0
    else
      vof_nod(i) = 1.0
    end if
  end do


  call Grid % Save_Debug_Vtu(append="vof_nod",       &
                             scalar_node=vof_nod,    &
                             scalar_name="vof_nod")

  !-----------------------!
  !   Find vof in cells   !
  !-----------------------!
  do c = 1, Grid % n_cells
    vof_cel(c) = 0.0
    do i_nod = 1, Grid % cells_n_nodes(c)
      i = Grid % cells_n(i_nod, c)
      vof_cel(c) = vof_cel(c) + vof_nod(i)
    end do
    vof_cel(c) = vof_cel(c) / real(Grid % cells_n_nodes(c))
  end do
  ! Exchange is needed here
  call Grid % Save_Debug_Vtu(append="vof_cel",       &
                             scalar_cell=vof_cel,    &
                             scalar_name="vof_cel")

  !------------------------------!
  !   Browse through cells ...   !
  !------------------------------!
  do c = 1, Grid % n_cells

    cut_count = 0

    call Polyhedron % Extract_From_Grid(Grid, c)
    Polyhedron % phi_int(1:Polyhedron % n_nodes) = UNDEFINED

    !----------------------------------------!
    !   Then browse through STL facets ...   !
    !----------------------------------------!
    do fac = 1, Sphere % n_facets

      ! STL vertex coordinates
      p1(1)=Sphere % x(1,fac); p1(2)=Sphere % y(1,fac); p1(3)=Sphere % z(1,fac)
      p2(1)=Sphere % x(2,fac); p2(2)=Sphere % y(2,fac); p2(3)=Sphere % z(2,fac)
      p3(1)=Sphere % x(3,fac); p3(2)=Sphere % y(3,fac); p3(3)=Sphere % z(3,fac)

      ! STL facet centroid
      f = (p1+p2+p3)/3.0

      ! STL facet normal
      n(1)=Sphere % nx(fac);  n(2)=Sphere % ny(fac);  n(3)=Sphere % nz(fac)

      ij_cut(:,:) = 0

      !----------------------------------------!
      !   Then browse through STL facets ...   !
      !----------------------------------------!
      do s = 1, Polyhedron % n_faces

        new_faces_n_nodes = 0  ! initalize number of nodes in new face
        new_faces_n(:)    = 0  ! initialize new face's node list to zero

        ! Browse nodes in circular direction
        do i_nod = 1, Polyhedron % faces_n_nodes(s)
          j_nod = i_nod + 1; if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1

          i = Polyhedron % faces_n(s, i_nod)
          j = Polyhedron % faces_n(s, j_nod)

          ! Plain copy of the old node (i) to new face's nodes
          new_faces_n_nodes = new_faces_n_nodes + 1
          new_faces_n(new_faces_n_nodes) = i  ! just copy "i"

          !---------------------------------------------------!
          !   The connection between i and j wasn't cut yet   !
          !---------------------------------------------------!
          if(ij_cut(i,j) == 0) then

            qi(1:3) = Polyhedron % nodes_xyz(i, 1:3)
            qj(1:3) = Polyhedron % nodes_xyz(j, 1:3)

            ! Do i and j cross STL facet?
            ij_cut_flag = .false.
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
                  ij_cut_flag = .true.
                end if
              end if  ! vol_3 and vol_4 have the same sign
            end if  ! vol_1 and vol_2 have different signs

            ! Yes, i and j cross the face !
            if(ij_cut_flag) then

              ! Vector connecting nodes i and j
              l = qj - qi

              ! Distance from i to a point on the plane
              ! (I picked center of the facet)
              di = dot_product(f-qi, n) / dot_product(l, n)

              cut_count = cut_count + 1

              ! Intersection point
              Polyhedron % n_nodes = Polyhedron % n_nodes + 1
              Polyhedron % nodes_xyz(Polyhedron % n_nodes, 1:3)  &
                                      = qi(1:3) + di * l(1:3)
              Polyhedron % phi_int(Polyhedron % n_nodes) = MEDIUM
              ! print '(a,i2,a,a,i2,a,i2)',                                     &
              !          ' # Inserting new node (', Polyhedron % n_nodes, ')',  &
              !          ' between nodes: ', i, ' and ', j

              ! Store new node which is added to intersection
              ij_cut(i,j) = Polyhedron % n_nodes
              ij_cut(j,i) = Polyhedron % n_nodes

              ! Add this new bloody node to the list of nodes in the face
              new_faces_n_nodes = new_faces_n_nodes + 1
              new_faces_n(new_faces_n_nodes) = ij_cut(i,j)

              ! Mark nodes which are on either side of interface ...
              ! ... using the scalar product with STL facet centre
              if( dot_product(n, qi-f) > 0.0 ) then
                Polyhedron % phi_int(i) = HIGH
              else
                Polyhedron % phi_int(i) = LOW
              end if
              if( dot_product(n, qj-f) > 0.0 ) then
                Polyhedron % phi_int(j) = HIGH
              else
                Polyhedron % phi_int(j) = LOW
              end if
            end if   ! ij_cut_flag

          !---------------------------------------------------!
          !   The connection between i and j was cut before   !
          !---------------------------------------------------!
          else  ! ij_cut(i,j) .ne 0

            ! print '(a,i2,a,a,i2,a,i2)',                        &
            !          ' #     The node (', ij_cut(i,j), ')',    &
            !          ' iss already inserted between nodes: ',  &
            !          i, ' and ', j

            ! But still, you have to add it to the face list
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n(new_faces_n_nodes) = ij_cut(i,j)

          end if  ! ij_cut == 0

        end do  ! through nodes of the face

        ! Now when all the nodes have been browsed, reform the
        ! face, I mean overwrite the old one with the new one
        Polyhedron % faces_n_nodes(s) = new_faces_n_nodes
        Polyhedron % faces_n(s, 1:new_faces_n_nodes)  &
                  = new_faces_n(1:new_faces_n_nodes)

      end do    ! through polyhedron faces
    end do      ! through STL facets

    !----------------------------------------------------------!
    !   Color the remaining polyhedron nodes by a flood fill   !
    !----------------------------------------------------------!
    if(cut_count .gt. 2) then
      do run = 1, Polyhedron % n_faces

        has_one_point_five = .false.

        do s = 1, Polyhedron % n_faces
          do i_nod = 1, Polyhedron % faces_n_nodes(s)
            j_nod = i_nod + 1; if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1

            i = Polyhedron % faces_n(s, i_nod)
            j = Polyhedron % faces_n(s, j_nod)

            ! Try not to spread across interface - if either i or j are
            ! newly formed nodes, their neighbours are alread set.
            if( Polyhedron % phi_int(i) .ne. MEDIUM .and.  &
                Polyhedron % phi_int(j) .ne. MEDIUM ) then

              if( Polyhedron % phi_int(i) .eq. UNDEFINED ) then
                has_one_point_five = .true.
                if( Polyhedron % phi_int(j) .eq. LOW ) then
                  Polyhedron % phi_int(i) = LOW
                end if
                if( Polyhedron % phi_int(j) .eq. HIGH ) then
                  Polyhedron % phi_int(i) = HIGH
                end if
              end if

              if( Polyhedron % phi_int(j) .eq. UNDEFINED ) then
                has_one_point_five = .true.
                if( Polyhedron % phi_int(i) .eq. LOW ) then
                  Polyhedron % phi_int(j) = LOW
                end if
                if( Polyhedron % phi_int(i) .eq. HIGH ) then
                  Polyhedron % phi_int(j) = HIGH
                end if
              end if

            end if

          end do  ! i, j, nodes
        end do    ! s, faces of polyhedron

        ! Exit if there is nothing left to color
        if(.not. has_one_point_five) then
          goto 1
        end if

      end do  ! run
    end if
1   continue

    !----------------------------------------------------------------!
    !                                                                !
    !   This is tough ... try to skip nodes with value 0.0 and 1.0   !
    !                                                                !
    !----------------------------------------------------------------!
    if(cut_count .gt. 2) then

      !-------------------------------------!
      !   First just copy the polyhedrons   !
      !-------------------------------------!
      call Pol0 % Create_From_Polyhedron(Polyhedron)
      call Pol1 % Create_From_Polyhedron(Polyhedron)

      !-----------------!
      !   Create Pol0   !
      !-----------------!
      do s = 1, Pol0 % n_faces
        new_faces_n_nodes = 0  ! initalize number of nodes in new face
        new_faces_n(:)    = 0  ! initialize new face's node list to zero

        ! Browse nodes in circular direction
        do i_nod = 1, Pol0 % faces_n_nodes(s)
          i = Pol0 % faces_n(s, i_nod)
          if( Pol0 % phi_int(i) .le. MEDIUM ) then
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n  (new_faces_n_nodes) = i  ! just copy "i"
          end if
        end do  ! through nodes of the face

        ! Now when all the nodes have been browsed, reform the
        ! face, I mean overwrite the old one with the new one
        Pol0 % faces_n_nodes(s) = new_faces_n_nodes
        Pol0 % faces_n(s,1:new_faces_n_nodes)  &
           = new_faces_n(1:new_faces_n_nodes)
      end do    ! through polyhedron faces
      ! What if there is a face without any nodes at MEDIUM?
      ! Well, faces_n_nodes will be zero for that face

      !-----------------------------------!
      !   Try to re-create missing face   !
      !-----------------------------------!
      n_combo     = 0
      combos(:,:) = 0
      do s = 1, Pol0 % n_faces
        do i_nod = 1, Pol0 % faces_n_nodes(s)
          j_nod = i_nod + 1; if(j_nod > Pol0 % faces_n_nodes(s)) j_nod = 1

          i = Pol0 % faces_n(s, i_nod)
          j = Pol0 % faces_n(s, j_nod)

          if( Pol0 % phi_int(i) .eq. MEDIUM .and.  &
              Pol0 % phi_int(j) .eq. MEDIUM ) then
            n_combo = n_combo + 1
            combos(n_combo, 1) = i
            combos(n_combo, 2) = j
            ! PRINT *, 'COMBO ', COMBOS(n_COMBO, :)
          end if
        end do
      end do

      if(n_combo > 0) then

        ! First two nodes are just coppied
        new_faces_n(1) = combos(1, 1)
        new_faces_n(2) = combos(1, 2)
        new_faces_n_nodes = 2
        next_n = new_faces_n(new_faces_n_nodes)
4       continue
        do cc = 2, n_combo
          if(combos(cc, 1) .eq. next_n) then
            next_n = combos(cc, 2)
            if(next_n .eq. new_faces_n(1)) goto 5  ! you closed the circle
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n(new_faces_n_nodes) = next_n
            goto 4
          end if
        end do
5       continue

        ! Add the face now
        Pol0 % n_faces       = Pol0 % n_faces + 1
        Pol0 % faces_n_nodes(Pol0 % n_faces) = new_faces_n_nodes
        Pol0 % faces_n(Pol0 % n_faces, 1:new_faces_n_nodes)  &
                         = new_faces_n(1:new_faces_n_nodes)
      end if

      !-----------------!
      !   Create Pol1   !
      !-----------------!
      do s = 1, Pol1 % n_faces
        new_faces_n_nodes = 0  ! initalize number of nodes in new face
        new_faces_n(:)    = 0  ! initialize new face's node list to zero

        ! Browse nodes in circular direction
        do i_nod = 1, Pol1 % faces_n_nodes(s)
          i = Pol1 % faces_n(s, i_nod)
          if( Pol1 % phi_int(i) .ge. MEDIUM) then
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n  (new_faces_n_nodes) = i  ! just copy "i"
          end if
        end do  ! through nodes of the face

        ! Now when all the nodes have been browsed, reform the
        ! face, I mean overwrite the old one with the new one
        Pol1 % faces_n_nodes(s) = new_faces_n_nodes
        Pol1 % faces_n(s,1:new_faces_n_nodes)  &
           = new_faces_n(1:new_faces_n_nodes)
      end do    ! through polyhedron faces
      ! What if there is a face without any nodes at MEDIUM?
      ! Well, faces_n_nodes will be zero for that face

      !-----------------------------------!
      !   Try to re-create missing face   !
      !-----------------------------------!
      n_combo     = 0
      combos(:,:) = 0
      do s = 1, Pol1 % n_faces
        do i_nod = 1, Pol1 % faces_n_nodes(s)
          j_nod = i_nod + 1; if(j_nod > Pol1 % faces_n_nodes(s)) j_nod = 1

          i = Pol1 % faces_n(s, i_nod)
          j = Pol1 % faces_n(s, j_nod)

          if( Pol1 % phi_int(i) .eq. MEDIUM .and.  &
              Pol1 % phi_int(j) .eq. MEDIUM ) then
            n_combo = n_combo + 1
            combos(n_combo, 1) = i
            combos(n_combo, 2) = j
            ! PRINT *, 'COMBO ', COMBOS(n_COMBO, :)
          end if
        end do
      end do

      if(n_combo > 0) then

        ! First two nodes are just coppied
        new_faces_n(1) = combos(1, 1)
        new_faces_n(2) = combos(1, 2)
        new_faces_n_nodes = 2
        next_n = new_faces_n(new_faces_n_nodes)
6       continue
        do cc = 2, n_combo
          if(combos(cc, 1) .eq. next_n) then
            next_n = combos(cc, 2)
            if(next_n .eq. new_faces_n(1)) goto 7  ! you closed the circle
            new_faces_n_nodes = new_faces_n_nodes + 1
            new_faces_n(new_faces_n_nodes) = next_n
            goto 6
          end if
        end do
7       continue

        ! Add the face now
        Pol1 % n_faces       = Pol1 % n_faces + 1
        Pol1 % faces_n_nodes(Pol1 % n_faces) = new_faces_n_nodes
        Pol1 % faces_n(Pol1 % n_faces, 1:new_faces_n_nodes)  &
                         = new_faces_n(1:new_faces_n_nodes)
      end if

    end if  ! cut count

    if(cut_count .gt. 2) then
      do i = 1, Polyhedron % n_nodes
        Polyhedron % phi(i) = real(Polyhedron % phi_int(i)) / real(HIGH)
        Pol0       % phi(i) = real(Pol0       % phi_int(i)) / real(HIGH)
        Pol1       % phi(i) = real(Pol1       % phi_int(i)) / real(HIGH)
      end do
      print *, '# Saving cell ', c, ' with ', cut_count, ' cuts'
      ! call Polyhedron % Plot_Polyhedron_Vtk("geo", c)
      call Pol0 % Plot_Polyhedron_Vtk("pol0", c)
      call Pol1 % Plot_Polyhedron_Vtk("pol1", c)
    end if

  end do  ! through cells

  call Work % Disconnect_Real_Node(vof_nod)
  call Work % Disconnect_Real_Cell(vof_cel)

  ! Exit gracefully
  call Comm_Mod_End
  stop

  end subroutine
