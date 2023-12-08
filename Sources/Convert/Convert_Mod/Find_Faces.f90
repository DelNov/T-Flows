!==============================================================================!
  subroutine Find_Faces(Convert, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is dedicated to finding faces within a grid, specifically
!>  for conformal meshes composed of tetrahedra, wedges, prisms, and hexahedra.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initialization: Allocates memory for arrays to store face and cell       !
!     information. It then calculates the maximum number of faces (max_fac)    !
!     and nodes (max_nod) per cell.                                            !
!   * Face identification in cells: The subroutine iterates over each cell,    !
!     identifying potential faces based on the nodes that define each cell.    !
!     It stores the smallest three nodes of each potential face, along with    !
!     the corresponding cell.                                                  !
!   * Sorting and annotating faces: The potential faces are sorted based on    !
!     their node indices. Faces with the same nodes are annotated for further  !
!     processing.                                                              !
!   * Main loop for face organization: Iterates through the annotated faces,   !
!     pairing cells to form internal faces of the mesh. It checks for matching !
!     nodes between cells to identify shared faces.                            !
!   * Finalizing face data: stores the identified faces' data in the Grid      !
!     structure, including the number of nodes per face and the cells each     !
!     face connects.                                                           !
!   * Memory Deallocation: deallocate local arrays to clean up the memory.     !
!------------------------------------------------------------------------------!
!   Note 1: To be more specific, this subrotine determines:                    !
!   * Grid % n_faces       - final number of faces (boundary + inside)         !
!   * Grid % faces_n_nodes - number of nodes for each face                     !
!   * Grid % faces_n       - nodes of each face                                !
!   * Grid % faces_c       - pair of cells surrounding each face               !
!   Note 2: Boundary faces have been determined in "Grid_Topology".            !
!   Note 3: This is called before geometrical quantities are calculated.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Convert_Type) :: Convert  !! parent class
  type(Grid_Type)     :: Grid     !! grid being converted
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, n1, n2, n3, f_nod(4)
  integer              :: n_match, i_fac, max_fac, match_nodes(-1:8)
  integer              :: i1, i2, i_nod, max_nod, cnt
  integer              :: fn(6,4)
  integer, allocatable :: face_n1(:)
  integer, allocatable :: face_n2(:)
  integer, allocatable :: face_n3(:)
  integer, allocatable :: face_cell(:), starts(:), ends(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Convert)
!==============================================================================!

  call Profiler % Start('Find_Faces')

  allocate(face_n1(Grid % n_cells*6))
  allocate(face_n2(Grid % n_cells*6))
  allocate(face_n3(Grid % n_cells*6))
  face_n1(:) = HUGE_INT
  face_n2(:) = HUGE_INT
  face_n3(:) = HUGE_INT

  allocate(face_cell(Grid % n_cells*6));  face_cell(:) = 0
  allocate(starts   (Grid % n_cells*6));  starts   (:) = 0
  allocate(ends     (Grid % n_cells*6));  ends     (:) = 0

  !---------------------------------------------------!
  !   Find max number of cells' faces for this grid   !
  !---------------------------------------------------!
  max_fac = 0
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .eq. 4) max_fac = max(max_fac, 4)
    if(Grid % cells_n_nodes(c) .eq. 5) max_fac = max(max_fac, 5)
    if(Grid % cells_n_nodes(c) .eq. 6) max_fac = max(max_fac, 5)
    if(Grid % cells_n_nodes(c) .eq. 8) max_fac = max(max_fac, 6)
  end do
  Assert(max_fac .ne. 0)

  !---------------------------------------------------!
  !   Find max number of faces' nodes for this grid   !
  !---------------------------------------------------!
  max_nod = 0
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .eq. 4) max_nod = max(max_nod, 3)
    if(Grid % cells_n_nodes(c) .eq. 5) max_nod = max(max_nod, 4)
    if(Grid % cells_n_nodes(c) .eq. 6) max_nod = max(max_nod, 4)
    if(Grid % cells_n_nodes(c) .eq. 8) max_nod = max(max_nod, 4)
  end do
  Assert(max_nod .ne. 0)

  !---------------------------------------------------!
  !   Fill the generic coordinates with some values   !
  !---------------------------------------------------!
  do c = 1, Grid % n_cells
    if(Grid % cells_n_nodes(c) .eq. 4) fn = TET
    if(Grid % cells_n_nodes(c) .eq. 5) fn = PYR
    if(Grid % cells_n_nodes(c) .eq. 6) fn = WED
    if(Grid % cells_n_nodes(c) .eq. 8) fn = HEX
    do i_fac = 1, max_fac
      if(Grid % cells_bnd_region(i_fac, c) .eq. 0) then

        ! Fetch face nodes (-1 becomes HUGE_INT and it will be the biggest)
        do i_nod = 1, max_nod
          f_nod(i_nod) = HUGE_INT
          if(fn(i_fac, i_nod) > 0) then
            f_nod(i_nod) = Grid % cells_n(fn(i_fac, i_nod), c)
          end if
        end do

        ! Sort them ...
        call Sort % Int_Array(f_nod(1:4))

        ! ... and store three of the smallest
        face_n1((c-1)*6 + i_fac) = f_nod(1)
        face_n2((c-1)*6 + i_fac) = f_nod(2)
        face_n3((c-1)*6 + i_fac) = f_nod(3)

        ! Store the cell too
        face_cell((c-1)*6 + i_fac) = c
      end if
    end do
  end do

  !--------------------------------------------------!
  !   Sort the cell faces according to coordinates   !
  !--------------------------------------------------!
  call Sort % Three_Int_Carry_Int(face_n1(1:Grid % n_cells*6),    &
                                  face_n2(1:Grid % n_cells*6),    &
                                  face_n3(1:Grid % n_cells*6),    &
                                  face_cell(1:Grid % n_cells*6))

  !------------------------------------------------!
  !   Anotate cell faces with same coordinates     !
  !------------------------------------------------!
  cnt = 1
  starts(1) = 1
  do c = 2, Grid % n_cells * 6
    if( face_n1(c) .ne. face_n1(c-1) .and.  &
        face_n2(c) .ne. face_n2(c-1) .and.  &
        face_n3(c) .ne. face_n3(c-1) ) then
      cnt = cnt + 1
      starts(cnt) = c
      ends(cnt-1) = c-1
    end if
  end do

  !---------------------------------------------!
  !                                             !
  !   Main loop to fill the faces_c structure   !
  !                                             !
  !---------------------------------------------!
  call Enlarge % Matrix_Int(Grid % cells_c, i=(/1,max_fac/))
  call Enlarge % Matrix_Int(Grid % faces_n, i=(/1,max_nod/))
  Assert(maxval(maxval(Grid % cells_c, dim=2), dim=1) == 0)

  do n3 = 1, cnt
    if(starts(n3) .ne. ends(n3)) then
      do i1=starts(n3),ends(n3)
        do i2=i1+1,ends(n3)

          !------------------------------------------!
          !   This is where it is set that c1 < c2   !
          !------------------------------------------!
          c1 = min(face_cell(i1),face_cell(i2))
          c2 = max(face_cell(i1),face_cell(i2))
          if(c1 .ne. c2) then

            !------------------------------!
            !   Number of matching nodes   !
            !------------------------------!
            n_match     = 0
            match_nodes = 0
            do n1 = 1, Grid % cells_n_nodes(c1)
              do n2 = 1, Grid % cells_n_nodes(c2)
                if(Grid % cells_n(n1,c1) .eq. Grid % cells_n(n2,c2)) then
                  n_match = n_match + 1
                  match_nodes(n1) = 1
                end if
              end do
            end do

            !-----------------------!
            !   general + general   !
            !     c1        c2      !
            !-----------------------!
            if(n_match > 2) then
              if(Grid % cells_n_nodes(c1) .eq. 4) fn = TET
              if(Grid % cells_n_nodes(c1) .eq. 5) fn = PYR
              if(Grid % cells_n_nodes(c1) .eq. 6) fn = WED
              if(Grid % cells_n_nodes(c1) .eq. 8) fn = HEX
              do i_fac = 1, max_fac
                if(Grid % cells_c(i_fac, c1) .eq. 0  .and.   & ! not set yet
                    ( max( match_nodes(fn(i_fac, 1)),0 ) + &
                      max( match_nodes(fn(i_fac, 2)),0 ) + &
                      max( match_nodes(fn(i_fac, 3)),0 ) + &
                      max( match_nodes(fn(i_fac, 4)),0 ) .eq. n_match ) ) then
                  Grid % n_faces = Grid % n_faces + 1
                  Grid % faces_c(1,Grid % n_faces) = c1
                  Grid % faces_c(2,Grid % n_faces) = c2
                  Grid % faces_n_nodes(Grid % n_faces) = n_match
                  do i_nod = 1, max_nod
                    if(fn(i_fac, i_nod) > 0) then
                      Grid % faces_n(i_nod,Grid % n_faces) =  &
                      Grid % cells_n(fn(i_fac, i_nod), c1)
                    end if
                  end do
                  Grid % cells_c(i_fac, c1) = 1  ! -> means: set
                end if
              end do
            end if   ! n_match .ne. 2
          end if   ! c1 .ne. c2
        end do   ! i2
      end do   ! i1
    end if
  end do    ! do n3

  print '(a38,i9)', '# Number of faces:                   ', Grid % n_faces

  ! De-allocate local memory
  deallocate(face_n1)
  deallocate(face_n2)
  deallocate(face_n3)
  deallocate(face_cell)
  deallocate(starts)
  deallocate(ends)

  call Profiler % Stop('Find_Faces')

  end subroutine
