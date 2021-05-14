!==============================================================================!
  subroutine Grid_Mod_Sort_Faces_Smart(grid)
!------------------------------------------------------------------------------!
!   Sorts array of faces in a smart way.  That would mean boundary faces       !
!   first, boundary region by boundary region, then inside faces, then         !
!   also according to indices of cells surrounding a face (c1 and c2).         !
!                                                                              !
!   Note that it doesn't store shadow faces (which are stored after normal     !
!   faces), but updates references to them.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(grid_type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, n, c, c1, c2, n_bc, color
  integer, allocatable :: old_nn  (:)
  integer, allocatable :: old_shad(:)
  integer, allocatable :: old_nods(:,:)
  real,    allocatable :: old_bxyz(:,:)
  integer, allocatable :: criteria(:,:)
  integer, parameter   :: BIG     = 2147483647  ! Euler's prime number
  logical, parameter   :: VERBOSE = .false.
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: M = MAX_FACES_N_NODES
!==============================================================================!

  ! Allocate memory
  allocate(criteria(grid % n_faces, 3))  ! 2nd ind smaller for memory alignemnt
  allocate(old_nn  (   grid % n_faces))  ! old number of nodes
  allocate(old_shad(   grid % n_faces))
  allocate(old_nods(M, grid % n_faces))
  allocate(old_bxyz(3,-grid % n_bnd_cells:-1))

  !--------------------------------------------------------!
  !   Form the three criteria:                             !
  !   1 - the strongest is boundary condition color; to    !
  !       sort faces by boundary condition colors first    !
  !   2 - the second on the list is cell number, so that   !
  !       inside cells are later browsed in a straight-    !
  !       forward incereasing way in the "Process"         !
  !   3 - third is the boundary cell number, as it will    !
  !       be completelly reformed soon anyway.             !
  !--------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    criteria(s,2) = c1
    criteria(s,3) = c2
    if(c2 > 0) then
      criteria(s,1) = BIG  ! make sure that inside faces end up last
    else
      criteria(s,1) = grid % bnd_cond % color(c2)
      criteria(s,3) = -criteria(s,3)
    end if
    grid % old_f(s) = s
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Int_Carry_Int(criteria(1:grid % n_faces, 1),  &
                                  criteria(1:grid % n_faces, 2),  &
                                  criteria(1:grid % n_faces, 3), grid % old_f)

  !------------------------------------------!
  !   Copy sorted c1s to face_c structure,   !
  !       and form c2 from the scratch       !
  !------------------------------------------!
  n_bc = 0
  do s = 1, grid % n_faces

    ! Copy c1 and c2 back ...
    grid % faces_c(1,s) = criteria(s,2)
    grid % faces_c(2,s) = criteria(s,3)

    ! ... but renumber c2 if on the boundary
    if(criteria(s,1) .ne. BIG) then         ! on the boundary
      n_bc = n_bc + 1                       ! increase the count
      grid % faces_c(2,s) = -n_bc           ! set face_c properly
      grid % old_c(-n_bc) = -criteria(s,3)  ! store the old number
    else
    end if
  end do

  !------------------------------------------------------!
  !     Using the old boundary cell number, retreive     !
  !   their boundary colors and geometrical quantities   !
  !------------------------------------------------------!
  do c=-1, -grid % n_bnd_cells, -1
    old_nn  ( -c) = grid % bnd_cond % color(grid % old_c(c))  ! use old_nn ...
    old_bxyz(1,c) = grid % xc(grid % old_c(c))                ! ... for colors
    old_bxyz(2,c) = grid % yc(grid % old_c(c))
    old_bxyz(3,c) = grid % zc(grid % old_c(c))
  end do
  do c=-1, -grid % n_bnd_cells, -1
    grid % bnd_cond % color(c) = old_nn  (-c)
    grid % xc(c) = old_bxyz(1,c)
    grid % yc(c) = old_bxyz(2,c)
    grid % zc(c) = old_bxyz(3,c)
  end do

  !---------------------------------------------!
  !   Sort faces_n_nodes, faces_n and faces_c   !
  !---------------------------------------------!
  do s = 1, grid % n_faces
    old_nn  (    s) = grid % faces_n_nodes( grid % old_f(s))
    old_nods(1:M,s) = grid % faces_n  (1:M, grid % old_f(s))
    old_shad(    s) = grid % faces_s      ( grid % old_f(s))
  end do
  do s = 1, grid % n_faces
    grid % faces_n_nodes(s) = old_nn  (     s)
    grid % faces_n(1:M,  s) = old_nods(1:M, s)
    grid % faces_s      (s) = old_shad(     s)
  end do

  !---------------------------------------------------!
  !   Copy topological quantities to boundary cells   !
  !---------------------------------------------------!
  do s = 1, grid % n_faces
    c2 = grid % faces_c(2, s)
    if(c2 < 0) then
      grid % cells_n_nodes(c2) = grid % faces_n_nodes(s)
      grid % cells_n(1:M,  c2) = grid % faces_n(1:M,  s)
    end if
  end do

  !------------------------------------------------!
  !   Sort geometriacal quantities for the faces   !
  !------------------------------------------------!

  ! Form the new face numbers from the old face numbers
  do s = 1, grid % n_faces
    grid % new_f(grid % old_f(s)) = s
  end do

  ! Do the actual sorting
  call Sort % Real_By_Index(grid % n_faces, grid % sx(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % sy(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % sz(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % dx(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % dy(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % dz(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % f (1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % xf(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % yf(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % zf(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % rx(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % ry(1), grid % new_f(1))
  call Sort % Real_By_Index(grid % n_faces, grid % rz(1), grid % new_f(1))

  ! Correct shadow faces
  do s = grid % n_faces + 1, grid % n_faces + grid % n_shadows
    grid % faces_s(s) = grid % new_f(grid % faces_s(s))
  end do

  ! Correct face indexes for cells
  do c = -grid % n_bnd_cells, grid % n_cells
    do n = 1, grid % cells_n_faces(c)
      s = grid % cells_f(n, c)             ! take the face's index
      grid % cells_f(n, c) = grid % new_f(s)
    end do
  end do

  ! Find boundary color ranges
  call Grid_Mod_Bnd_Cond_Ranges(grid)

  if(VERBOSE) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)
      if(c2 < 0) then
        print '(4i6)', s, c1, c2, grid % bnd_cond % color(c2)
      end if
    end do

    do color = 1, grid % n_bnd_cond
      print '(2i6)', grid % bnd_cond % color_s_cell(color),  &
                     grid % bnd_cond % color_e_cell(color)
    end do
  end if

  end subroutine
