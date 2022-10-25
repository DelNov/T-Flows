!==============================================================================!
  subroutine Sort_Faces_Smart(Grid)
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
  class(Grid_Type) :: Grid
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: BIG     = 2147483647  ! Euler's prime number
  logical, parameter :: VERBOSE = .false.
!-----------------------------------[Locals]-----------------------------------!
  integer              :: s, m, n, c, c1, c2, n_bc, color
  integer, allocatable :: old_nn  (:)
  integer, allocatable :: old_shad(:)
  integer, allocatable :: old_nods(:,:)
  real,    allocatable :: old_bxyz(:,:)
  integer, allocatable :: criteria(:,:)
!==============================================================================!

  ! Allocate memory
  m = size(Grid % faces_n, 1)
  allocate(criteria(Grid % n_faces, 3))  ! 2nd ind smaller for memory alignemnt
  allocate(old_nn  (   Grid % n_faces))  ! old number of nodes
  allocate(old_shad(   Grid % n_faces))
  allocate(old_nods(m, Grid % n_faces))
  allocate(old_bxyz(3,-Grid % n_bnd_cells:-1))

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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    criteria(s,2) = c1
    criteria(s,3) = c2
    if(c2 > 0) then
      criteria(s,1) = BIG  ! make sure that inside faces end up last
    else
      criteria(s,1) = Grid % bnd_cond % color(c2)
      criteria(s,3) = -criteria(s,3)
    end if
    Grid % old_f(s) = s
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort % Three_Int_Carry_Int(criteria(1:Grid % n_faces, 1),  &
                                  criteria(1:Grid % n_faces, 2),  &
                                  criteria(1:Grid % n_faces, 3), Grid % old_f)

  !------------------------------------------!
  !   Copy sorted c1s to face_c structure,   !
  !       and form c2 from the scratch       !
  !------------------------------------------!
  n_bc = 0
  do s = 1, Grid % n_faces

    ! Copy c1 and c2 back ...
    Grid % faces_c(1,s) = criteria(s,2)
    Grid % faces_c(2,s) = criteria(s,3)

    ! ... but renumber c2 if on the boundary
    if(criteria(s,1) .ne. BIG) then         ! on the boundary
      n_bc = n_bc + 1                       ! increase the count
      Grid % faces_c(2,s) = -n_bc           ! set face_c properly
      Grid % old_c(-n_bc) = -criteria(s,3)  ! store the old number
    else
    end if
  end do

  !------------------------------------------------------!
  !     Using the old boundary cell number, retreive     !
  !   their boundary colors and geometrical quantities   !
  !------------------------------------------------------!
  do c=-1, -Grid % n_bnd_cells, -1
    old_nn  ( -c) = Grid % bnd_cond % color(Grid % old_c(c))  ! use old_nn ...
    old_bxyz(1,c) = Grid % xc(Grid % old_c(c))                ! ... for colors
    old_bxyz(2,c) = Grid % yc(Grid % old_c(c))
    old_bxyz(3,c) = Grid % zc(Grid % old_c(c))
  end do
  do c=-1, -Grid % n_bnd_cells, -1
    Grid % bnd_cond % color(c) = old_nn  (-c)
    Grid % xc(c) = old_bxyz(1,c)
    Grid % yc(c) = old_bxyz(2,c)
    Grid % zc(c) = old_bxyz(3,c)
  end do

  !---------------------------------------------!
  !   Sort faces_n_nodes, faces_n and faces_c   !
  !---------------------------------------------!
  do s = 1, Grid % n_faces
    old_nn  (    s) = Grid % faces_n_nodes( Grid % old_f(s))
    old_nods(1:m,s) = Grid % faces_n  (1:m, Grid % old_f(s))
    old_shad(    s) = Grid % faces_s      ( Grid % old_f(s))
  end do
  do s = 1, Grid % n_faces
    Grid % faces_n_nodes(s) = old_nn  (     s)
    Grid % faces_n(1:m,  s) = old_nods(1:m, s)
    Grid % faces_s      (s) = old_shad(     s)
  end do

  !---------------------------------------------------!
  !   Copy topological quantities to boundary cells   !
  !---------------------------------------------------!
  do s = 1, Grid % n_faces
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      Grid % cells_n_nodes(c2) = Grid % faces_n_nodes(s)
      Grid % cells_n(1:Grid % cells_n_nodes(c2), c2) =  &
      Grid % faces_n(1:Grid % faces_n_nodes(s), s)
    end if
  end do

  !------------------------------------------------!
  !   Sort geometriacal quantities for the faces   !
  !------------------------------------------------!

  ! Form the new face numbers from the old face numbers
  do s = 1, Grid % n_faces
    Grid % new_f(Grid % old_f(s)) = s
  end do

  ! Do the actual sorting
  call Sort % Real_By_Index(Grid % n_faces, Grid % sx(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % sy(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % sz(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % dx(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % dy(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % dz(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % f (1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % xf(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % yf(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % zf(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % rx(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % ry(1), Grid % new_f(1))
  call Sort % Real_By_Index(Grid % n_faces, Grid % rz(1), Grid % new_f(1))

  ! Correct shadow faces
  do s = Grid % n_faces + 1, Grid % n_faces + Grid % n_shadows
    Grid % faces_s(s) = Grid % new_f(Grid % faces_s(s))
  end do

  ! Correct face indexes for cells
  do c = -Grid % n_bnd_cells, Grid % n_cells
    do n = 1, Grid % cells_n_faces(c)
      s = Grid % cells_f(n, c)             ! take the face's index
      Grid % cells_f(n, c) = Grid % new_f(s)
    end do
  end do

  ! Find boundary color ranges
  call Bnd_Cond_Ranges(Grid)

  if(VERBOSE) then
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      if(c2 < 0) then
        print '(4i6)', s, c1, c2, Grid % bnd_cond % color(c2)
      end if
    end do

    do color = 1, Grid % n_bnd_cond
      print '(2i6)', Grid % bnd_cond % color_s_cell(color),  &
                     Grid % bnd_cond % color_e_cell(color)
    end do
  end if

  end subroutine
