!==============================================================================!
  subroutine Grid_Mod_Coarsen(grid)
!------------------------------------------------------------------------------!
!   Coarsens the grid with METIS library.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, nc, nc1, nc2, s, sl, lp, i, lev, lev_parts
  integer              :: n_cells, n_faces, n_parts, arr_s, arr_e, val_1, val_2
  integer              :: c_lev
  integer, allocatable :: c1_arr(:), c2_arr(:)
  integer, allocatable :: cells_c(:,:), cells_n_cells(:), faces_c(:,:),  &
                          new_c(:), old_c(:),                            &
                          cell_mapping(:,:)

! Variabes for call to METIS
  integer              :: n_constrains,      &  ! number of constraints
                          return_val            ! return value from METIS
  integer, allocatable :: row(:),            &  ! incidency matrix in ...
                          col(:),            &  ! compresser row format
                          cell_weights(:),   &  ! weights of vertices
                          face_weights(:),   &  ! weights of edges
                          vert_data(:),      &  ! amount of data for vertices
                          part(:)               ! result of partitioning
  real, allocatable    :: part_weight(:),    &
                          imbalance(:)          ! allowed imbalance
  integer              :: metis_options(41)     ! options passed to METIS
!==============================================================================!

  call Grid_Mod_Allocate_Levels(grid)

  print *, '#==================================='
  print *, '# Coarsening the grid for multigrid '
  print *, '#-----------------------------------'

  ! Allocate memory
  allocate(new_c(grid % n_cells))
  allocate(old_c(grid % n_cells))

  ! Number of cells (vertices for METIS) and number of faces ...
  ! ... inside the grid (edges for METIS) for the first level
  n_cells = grid % n_cells
  n_faces = 0
  do s = 1, grid % n_faces
    c2 = grid % faces_c(2,s)
    if(grid % faces_c(2,s) > 0) n_faces = n_faces + 1
  end do

  ! Once n_cells and n_faces are known, allocate memory
  allocate(faces_c      ( 2,n_faces));  faces_c      (:,:) = 0
  allocate(cells_n_cells(   n_cells));  cells_n_cells(:)   = 0
  allocate(cells_c      (24,n_cells));  cells_c      (:,:) = 0
  allocate(part         (   n_cells));  part         (:)   = 1
  allocate(row          ( 1+n_cells))
  allocate(col          ( 2*n_faces))

  allocate(c1_arr(n_faces))
  allocate(c2_arr(n_faces))

  ! Allocate memory for control parameters for METIS
  n_constrains = 1
  n_parts      = 4   ! how many partitions at each level
  allocate(imbalance   (n_constrains))
  allocate(cell_weights(n_cells * n_constrains))
  allocate(face_weights(n_faces + n_faces))
  allocate(vert_data   (n_cells))
  allocate(part_weight (n_parts * n_constrains))

  !---------------------------------------------------------!
  !   Initialize arrays for level zero and coarser levels   !
  !---------------------------------------------------------!

  ! Cells on level 0
  grid % level(0) % n_cells = grid % n_cells
  do c = 1, grid % n_cells
    grid % level(0) % cell(c) = c
  end do

  ! Faces on level 0
  grid % level(0) % n_faces = grid % n_faces
  do s = 1, grid % n_faces
    grid % level(0) % face(s) = s
    grid % level(0) % faces_c(1, s) = grid % faces_c(1, s)
    grid % level(0) % faces_c(2, s) = grid % faces_c(2, s)
  end do

  ! Cells and faces on other levels
  do lev = 1, MAX_MG_LEV
    do c = 1, grid % n_cells
      grid % level(lev) % cell(c) = 1
    end do
    do s = 1, grid % n_faces
      grid % level(lev) % face(s) = 1
    end do
  end do

  ! Find out the number of effective levels
  do lev = 1, MAX_MG_LEV
    lev_parts = n_parts ** (lev-1)
    if(n_parts**lev > grid % n_cells/27) then
      grid % n_levels = lev
      goto 1
    end if
  end do
  grid % n_levels = MAX_MG_LEV
1 continue

  print '(a29,i2,a8)', ' # Will perform coarsening in',  &
                       grid % n_levels, ' levels.'

  !---------------------------!
  !                           !
  !   Browse through levels   !
  !                           !
  !---------------------------!
  do lev = 1, grid % n_levels

    lev_parts = n_parts ** (lev-1)

    ! Where to store this level (sl)
    sl = grid % n_levels - lev + 1
    print '(a19,i2)', ' # Working on level', sl

    !---------------------------------------------!
    !   Through exising partitions on level lev   !
    !---------------------------------------------!
    new_c(:) = 0
    old_c(:) = 0
    do lp = 1, lev_parts

      ! Mark cells in current partition
      i = 0
      do c = 1, grid % n_cells
        if(grid % level(sl+1) % cell(c) .eq. lp) then
          i = i + 1
          new_c(c) = i
          old_c(i) = c
        end if
      end do
      n_cells = i

      ! Form face (edge for METIS) connectivity in ...
      ! ... current partition at current level
      faces_c(:,:) = 0
      i = 0
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1, s)
        c2 = grid % faces_c(2, s)

        if(c2 > 0) then
          if(grid % level(sl+1) % cell(c1) .eq. lp .and.  &
             grid % level(sl+1) % cell(c2) .eq. lp) then
            i = i + 1
            faces_c(1,i) = new_c(c1)
            faces_c(2,i) = new_c(c2)
          end if
        end if
      end do
      n_faces = i

      ! Form cell to cell connectivity for partition and level
      cells_n_cells(:) = 0
      cells_c   (:,:) = 0
      do s = 1, n_faces
        nc1 = faces_c(1,s)
        nc2 = faces_c(2,s)

        cells_n_cells(nc1) = cells_n_cells(nc1) + 1
        cells_n_cells(nc2) = cells_n_cells(nc2) + 1

        cells_c(cells_n_cells(nc1), nc1) = nc2
        cells_c(cells_n_cells(nc2), nc2) = nc1
      end do

      ! Fill-up the structures needed to call METIS (row and col)
      row(:) = 0
      row(1) = 0
      do nc = 1, n_cells
        row(nc+1) = row(nc) + cells_n_cells(nc)
      end do

      col(:) = 0
      do nc = 1, n_cells
        do i = 1, cells_n_cells(nc)
          col(row(nc) + i) = cells_c(i, nc) - 1  ! -1, METIS works from 0
        end do
      end do

      imbalance(:)    = 1.001
      cell_weights(:) = 1
      face_weights(:) = 1
      vert_data(:)    = 1
      part_weight(:)  = 1.0 / real(n_parts)

      !-----------------------------------------!
      !   Exectute the call to METIS function   !
      !-----------------------------------------!

      metis_options = -1                       ! Initialize all to default
      metis_options(METIS_OPTION_DBGLVL) = 0
      metis_options(METIS_OPTION_CONTIG) = 1

      call METIS_PartGraphRecursive(n_cells,       &  !  1. (in), int
                                    n_constrains,  &  !  2. (in), int
                                    row,           &  !  3. (in), int(:)
                                    col,           &  !  4. (in), int(:)
                                    cell_weights,  &  !  5. (in), int(:)
                                    vert_data,     &  !  6. (in), int(:)
                                    face_weights,  &  !  7. (in), int(:)
                                    n_parts,       &  !  8. (in), int(:)
                                    part_weight,   &  !  9. (in), real(:)
                                    imbalance,     &  ! 10. (in), real(:)
                                    metis_options, &  ! 11. (in), int(:)
                                    return_val,    &  ! 12. (out) int(:)
                                    part(:))          ! 13. (out) int(:)

      part(1:n_cells) = part(1:n_cells) + 1  ! +1, METIS works from zero

      !-----------------------------------------------------!
      !   Save the result from the call to METIS function   !
      !-----------------------------------------------------!
      if(lev .eq. 1) then
        do c = 1, n_cells
          grid % level(sl) % cell(c) = part(c)
        end do
      else
        do nc = 1, n_cells
          c = old_c(nc)
          grid % level(sl) % cell(c) =   &
           (grid % level(sl+1) % cell(c)-1) * n_parts + part(nc)
        end do
      end if
    end do
  end do

  !-------------------------------------------------!
  !                                                 !
  !   Do the necessary book-keeping on all levels   !
  !                                                 !
  !-------------------------------------------------!

  !---------------------------------------------!
  !     Find cells around faces and and ...     !
  !   ... count number of faces at each level   !
  !---------------------------------------------!

  do lev = 0, grid % n_levels
    i = 0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)
 
      if(c2 > 0) then
        c1 = grid % level(lev) % cell(c1)
        c2 = grid % level(lev) % cell(c2)

        if(c1 .ne. c2) then
          i = i + 1
          grid % level(lev) % faces_c(1, i) = min(c1, c2)
          grid % level(lev) % faces_c(2, i) = max(c1, c2)
        end if
      end if
    end do
    grid % level(lev) % n_faces = i
    print *, '# Non-compressed number of faces at level ',  &
             lev, '=', grid % level(lev) % n_faces
  end do

  !-------------------------------!
  !   Compress face information   !
  !-------------------------------!

  do lev = 0, grid % n_levels

    ! Store face information into spare arrays
    c1_arr(:) = 0
    c2_arr(:) = 0
    do s = 1, grid % level(lev) % n_faces
      c1_arr(s) = grid % level(lev) % faces_c(1, s)
      c2_arr(s) = grid % level(lev) % faces_c(2, s)
    end do

    ! Sort faces by both indexes
    call Sort_Mod_2_Int(c1_arr(1:grid % level(lev) % n_faces),  &
                        c2_arr(1:grid % level(lev) % n_faces))

    write(lev+200,*), 'Level: ', lev, 'after sorting faces by second index'
    do s = 1, grid % level(lev) % n_faces
      write(lev+200,'(3i6)'), s, c1_arr(s), c2_arr(s)
    end do

    ! Selece faces by both cell indices
    i     = 1
    val_1 = c1_arr(1)
    val_2 = c2_arr(1)
    grid % level(lev) % faces_c(1, i) = val_1
    grid % level(lev) % faces_c(2, i) = val_2

    do s = 2, grid % level(lev) % n_faces
      if(c1_arr(s) .ne. val_1 .or.  &
         c2_arr(s) .ne. val_2) then
        i = i + 1
        grid % level(lev) % faces_c(1, i) = c1_arr(s)
        grid % level(lev) % faces_c(2, i) = c2_arr(s)
        val_1 = c1_arr(s)
        val_2 = c2_arr(s)
      end if
    end do
    grid % level(lev) % n_faces = i
    print *, '# Compressed number of faces at level ',  &
             lev, '=', grid % level(lev) % n_faces

    write(lev+300,*), 'Level: ', lev, 'after compressing faces'
    do s = 1, grid % level(lev) % n_faces
      write(lev+300,'(3i6)'), s,                                  &
                              grid % level(lev) % faces_c(1, s),  &
                              grid % level(lev) % faces_c(2, s)
    end do

  end do

  !-----------------------------------------!
  !   Count number of cells at each level   !
  !-----------------------------------------!
  grid % level(0) % n_cells = grid % n_cells  ! level 0 is the non-coarsen one
  do lev = 1, grid % n_levels
    grid % level(lev) % n_cells = maxval(grid % level(lev) % cell(:))
  end do

  !-------------------------------------!
  !                                     !
  !   Check sanity of the grid levels   !
  !                                     !
  !-------------------------------------!
  i = 0
  do lev = 1, grid % n_levels
    i = max(i, grid % level(lev) % n_cells)
  end do
  allocate(cell_mapping(MAX_MG_LEV, i)); cell_mapping = 0

  do lev = 1, grid % n_levels - 1
    n_parts = maxval(grid % level(lev) % cell(:))
    print '(a,i2,a,i2)', ' # Checking levels', lev, ' and', lev+1

    ! Browse through parts of this level
    do lp = 1, n_parts
      do c = 1, grid % n_cells
        if(grid % level(lev) % cell(c) == lp) then
          if(cell_mapping(lev, lp) .eq. 0) then
            cell_mapping(lev, lp) = grid % level(lev+1) % cell(c)
          else
            if(cell_mapping(lev, lp) .ne. grid % level(lev+1) % cell(c)) then
              print *, '# Mapping failed at level ', lev
              print *, '# Stopping the program!   '
              exit
            end if
          end if
        end if
      end do
    end do

  end do

  !------------------------------------------------!
  !   Compute cell coordinates at coarser levels   !
  !------------------------------------------------!
  do lev = 0, grid % n_levels

    print *, ' level n_cells = ', grid % level(lev) % n_cells

    grid % level(lev) % xc(:) = 0
    grid % level(lev) % yc(:) = 0
    grid % level(lev) % zc(:) = 0
    grid % level(lev) % n_finest_cells(:) = 0

    do c = 1, grid % n_cells

      c_lev = grid % level(lev) % cell(c)
      grid % level(lev) % xc(c_lev) = grid % level(lev) % xc(c_lev)  &
                                    + grid % xc(c)
      grid % level(lev) % yc(c_lev) = grid % level(lev) % yc(c_lev)  &
                                    + grid % yc(c)
      grid % level(lev) % zc(c_lev) = grid % level(lev) % zc(c_lev)  &
                                    + grid % zc(c)

      grid % level(lev) % n_finest_cells(c_lev) =    &
      grid % level(lev) % n_finest_cells(c_lev) + 1
    end do

    do c_lev = 1, grid % level(lev) % n_cells
      grid % level(lev) % xc(c_lev) = grid % level(lev) % xc(c_lev)  &
                                    / grid % level(lev) % n_finest_cells(c_lev)
      grid % level(lev) % yc(c_lev) = grid % level(lev) % yc(c_lev)  &
                                    / grid % level(lev) % n_finest_cells(c_lev)
      grid % level(lev) % zc(c_lev) = grid % level(lev) % zc(c_lev)  &
                                    / grid % level(lev) % n_finest_cells(c_lev)
    end do
  end do

  deallocate(col)
  deallocate(row)
  deallocate(part_weight)
  deallocate(vert_data)
  deallocate(face_weights)
  deallocate(cell_weights)
  deallocate(imbalance)
  deallocate(old_c)
  deallocate(new_c)

  end subroutine
