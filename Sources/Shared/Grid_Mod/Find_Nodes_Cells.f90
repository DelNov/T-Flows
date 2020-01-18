!==============================================================================!
  subroutine Grid_Mod_Find_Nodes_Cells(grid)
!------------------------------------------------------------------------------!
!   Cells around each node are needed for Lagrangian particle tracking.        !
!                                                                              !
!   The suboutine browses through faces to determine node to cell  connect-    !
!   ivity.  Browsing through faces is needed in order to look for periodic     !
!   cells.  For cases which are periodic in two directions, faces have to      !
!   browse through first and second neighbours.  For that, we need cell neig-  !
!   hbours.                                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, lc1, lc2, c1, c2, ln1, ln2, n1, n2, s, n, ps
  integer              :: run
  integer              :: max_n_cells       ! max number of cells around nodes
  real                 :: lx, ly, lz, x1, y1, z1, x2, y2, z2
  logical              :: px, py, pz
  integer, allocatable :: nodes_c(:,:)      ! local storage of nodes' cells
  integer, allocatable :: cells_n_cells(:)  ! local storage of cells neighbrs
  integer, allocatable :: cells_c    (:,:)
  integer, allocatable :: faces_n_cells(:)  ! local storage for faces' cells
  integer, allocatable :: faces_c(:,:)      ! can hold 1st or 2nd neighbours
  integer              :: work(128)
!==============================================================================!

  ! Determine cells neighbours
  allocate(cells_n_cells(-grid % n_bnd_cells:grid % n_cells))
  allocate(faces_n_cells( grid % n_faces));  faces_n_cells(:) = 0

  !------------------------------------!
  !   Find cells' neighbouring cells   !
  !------------------------------------!

  do run = 1, 2
    cells_n_cells(:) = 0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      ! Count maximum number of cells' neighbours
      cells_n_cells(c1) = cells_n_cells(c1) + 1
      cells_n_cells(c2) = cells_n_cells(c2) + 1
      ! Store cells' neighbours
      if(run .eq. 2) then
        cells_c(cells_n_cells(c1), c1) = c2
        cells_c(cells_n_cells(c2), c2) = c1
      end if
    end do

    ! Allocate memory for cells' neighbours
    if(run .eq. 1) then
      max_n_cells = maxval(cells_n_cells(:))
      allocate(cells_c(max_n_cells, -grid % n_bnd_cells:grid % n_cells))
      cells_c(:,:) = 0
    end if
  end do

  !-------------------------------------------------------------------!
  !   Find cells which surround faces locally (for this subroutine)   !
  !-------------------------------------------------------------------!

  ! Allocate memory
  allocate(faces_c(2*max_n_cells, grid % n_faces))

  ! For all cells, take only first cells
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    faces_n_cells(s) = 2
    faces_c(1, s) = c1
    faces_c(2, s) = c2
  end do

  ! Extend to 2nd neighbours for periodic faces
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if( Math_Mod_Distance(grid % xc(c1), grid % yc(c1), grid % zc(c1),    &
                          grid % xc(c2), grid % yc(c2), grid % zc(c2))    &
        >                                                                 &
        2.0 * sqrt(grid % dx(s)**2 + grid % dy(s)**2 + grid % dz(s)**2) ) then
      work(1:2*max_n_cells) = HUGE_INT
      n1 = cells_n_cells(c1)
      n2 = cells_n_cells(c2)
      work(   1:n1   ) = cells_c(1:n1,c1)
      work(n1+1:n1+n2) = cells_c(1:n2,c2)
      n = n1 + n2
      call Sort_Mod_Unique_Int(n, work(1:n))
      faces_n_cells(s) = n
      faces_c (1:n, s) = work(1:n)
    end if
  end do
  deallocate(cells_c)
  deallocate(cells_n_cells)

  ! Allocate memory for node cells
  n = grid % n_nodes
  allocate(grid % nodes_n_cells(1:n))
  grid % nodes_n_cells(:) = 0

  ! Take aliases
  px = .false.
  py = .false.
  pz = .false.
  lx = grid % per_x;  if(abs(lx) > TINY) px = .true.
  ly = grid % per_y;  if(abs(ly) > TINY) py = .true.
  lz = grid % per_z;  if(abs(lz) > TINY) pz = .true.

  ! This information is missing :-(
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      grid % cells_n_nodes(c2) = grid % faces_n_nodes(s)
      do n = 1, grid % cells_n_nodes(c2)
        grid % cells_n(n,c2) = grid % faces_n(n,s)
      end do
    end if
  end do

  !--------------------------------------------------------!
  !   Find maximum number of cells surrounding each node   !
  !--------------------------------------------------------!

  do run = 1, 2

  ! (Re)initialize the cell count
  grid % nodes_n_cells(:) = 0

  ! All faces
  do s = 1, grid % n_faces
    do lc1 = 1, faces_n_cells(s)
      do lc2 = lc1+1, faces_n_cells(s)
        c1 = faces_c(lc1,s)
        c2 = faces_c(lc2,s)
        do ln1 = 1, grid % cells_n_nodes(c1)    ! local node number 1
          n1 = grid % cells_n(ln1, c1)          ! global node number 1
          x1 = grid % xn(n1)
          y1 = grid % yn(n1)
          z1 = grid % zn(n1)
          do ln2 = 1, grid % cells_n_nodes(c2)  ! local node number 2
            n2 = grid % cells_n(ln2, c2)        ! global node number 2
            x2 = grid % xn(n2)
            y2 = grid % yn(n2)
            z2 = grid % zn(n2)
            if(Math_Mod_Distance_Squared(x1, y1, z1, x2, y2, z2) < PICO) then
              if(run .eq. 1) then
                grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
                grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
              end if
              if(run .eq. 2) then
                grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                nodes_c(grid % nodes_n_cells(n1), n1) = c1
                grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                nodes_c(grid % nodes_n_cells(n1), n1) = c2
                grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                nodes_c(grid % nodes_n_cells(n2), n2) = c1
                grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                nodes_c(grid % nodes_n_cells(n2), n2) = c2
              end if
            end if
            if( px .neqv. py .neqv. pz ) then  ! only one periodic direction
              if( ( px .and.                                             &
                    (Math_Mod_Distance_Squared(x1,   y1,z1,              &
                                               x2+lx,y2,z2) < PICO .or.  &
                     Math_Mod_Distance_Squared(x1,   y1,z1,              &
                                               x2-lx,y2,z2) < PICO)      &
                  ) .or.                                                 &
                  ( py .and.                                             &
                    (Math_Mod_Distance_Squared(x1,y1,   z1,              &
                                               x2,y2+ly,z2) < PICO .or.  &
                     Math_Mod_Distance_Squared(x1,y1,   z1,              &
                                               x2,y2-ly,z2) < PICO)      &
                  ) .or.                                                 &
                  ( pz .and.                                             &
                    (Math_Mod_Distance_Squared(x1,y1,z1,                 &
                                               x2,y2,z2+lz) < PICO .or.  &
                     Math_Mod_Distance_Squared(x1,y1,z1,                 &
                                               x2,y2,z2-lz) < PICO)      &
                  ) ) then
                if(run .eq. 1) then
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
                end if
                if(run .eq. 2) then
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                  nodes_c(grid % nodes_n_cells(n1), n1) = c1
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                  nodes_c(grid % nodes_n_cells(n1), n1) = c2
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                  nodes_c(grid % nodes_n_cells(n2), n2) = c1
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                  nodes_c(grid % nodes_n_cells(n2), n2) = c2
                end if
              end if
            end if
            if( px .and. py .and. .not. pz ) then  ! x and y are periodic
              if(Math_Mod_Distance_Squared(x1,   y1,z1,                 &
                                           x2+lx,y2,z2) < PICO .or.     &
                 Math_Mod_Distance_Squared(x1,   y1,z1,                 &
                                           x2-lx,y2,z2) < PICO .or.     &
                 Math_Mod_Distance_Squared(x1,y1,z1,                    &
                                           x2,y2+ly,z2) < PICO .or.     &
                 Math_Mod_Distance_Squared(x1,y1,z1,                    &
                                           x2,y2-ly,z2) < PICO .or.     &
                 Math_Mod_Distance_Squared(x1,   y1,   z1,              &
                                           x2+lx,y2+ly,z2) < PICO .or.  &
                 Math_Mod_Distance_Squared(x1,   y1,   z1,              &
                                           x2+lx,y2-ly,z2) < PICO .or.  &
                 Math_Mod_Distance_Squared(x1,   y1,   z1,              &
                                           x2-lx,y2+ly,z2) < PICO .or.  &
                 Math_Mod_Distance_Squared(x1,   y1,   z1,              &
                                           x2-lx,y2-ly,z2) < PICO) then
                if(run .eq. 1) then
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 2
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 2
                end if
                if(run .eq. 2) then
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                  nodes_c(grid % nodes_n_cells(n1), n1) = c1
                  grid % nodes_n_cells(n1) = grid % nodes_n_cells(n1) + 1
                  nodes_c(grid % nodes_n_cells(n1), n1) = c2
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                  nodes_c(grid % nodes_n_cells(n2), n2) = c1
                  grid % nodes_n_cells(n2) = grid % nodes_n_cells(n2) + 1
                  nodes_c(grid % nodes_n_cells(n2), n2) = c2
                end if
              end if
            end if
          end do
        end do
      end do
    end do
  end do

  ! Allocate memory for cells surrounding each node
  if(run .eq. 1) then
    max_n_cells = maxval(grid % nodes_n_cells)
    allocate(nodes_c(1:max_n_cells, 1:grid % n_nodes))
  end if

  end do

  deallocate(faces_n_cells)  ! local faces_n_cells not used beyond this point
  deallocate(faces_c)        ! local faces_c not used beyond this point

  do n = 1, grid % n_nodes
    call Sort_Mod_Unique_Int( grid % nodes_n_cells(n),  &
                nodes_c(1:grid % nodes_n_cells(n), n) )
  end do

  max_n_cells = maxval(grid % nodes_n_cells)

  ! Allocate memory for cells surrounding each node
  allocate(grid % nodes_c(1:max_n_cells, 1:grid % n_nodes))
  grid % nodes_c(:,:) = 0

  do n = 1, grid % n_nodes
    grid % nodes_c(1:grid % nodes_n_cells(n), n) = &
       nodes_c(1:grid % nodes_n_cells(n), n)
  end do

  deallocate(nodes_c)

  !----------------------------------------------!
  !   Store boundary cells to face connections   !
  !----------------------------------------------!
  do s = 1, grid % n_bnd_cells
    c2 = grid % faces_c(2,s)
    grid % cells_bnd_face(c2) = s
  end do

  end subroutine
