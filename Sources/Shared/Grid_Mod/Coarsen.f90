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
  integer              :: n_verts, n_edges, n_parts, arr_s, arr_e, val_1, val_2
  integer, allocatable :: arr_1(:), arr_2(:)
  integer, allocatable :: star(:,:), star_size(:),   &
                          edges_v(:,:), edges_c(:),  &
                          new_c(:), old_c(:),        &
                          cell_mapping(:,:)

! Variabes for call to METIS
  integer              :: n_constrains,      &  ! number of constraints
                          return_val            ! return value from METIS
  integer, allocatable :: row(:),            &  ! incidency matrix in ...
                          col(:),            &  ! compresser row format
                          vert_weights(:),   &  ! weights of vertices
                          edge_weights(:),   &  ! weights of edges
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

  ! Number of vertices and number of edges for first level 
  n_verts = grid % n_cells
  n_edges = 0
  do s = 1, grid % n_faces
    c2 = grid % faces_c(2,s)
    if(grid % faces_c(2,s) > 0) n_edges = n_edges + 1
  end do

  ! Once n_verts(1) and n_edegs(1) are known, allocate memory
  allocate(edges_v  ( 2,n_edges));  edges_v  (:,:) = 0
  allocate(edges_c  (   n_edges));  edges_c  (:)   = 0
  allocate(star_size(   n_verts));  star_size(:)   = 0
  allocate(star     (24,n_verts));  star     (:,:) = 0
  allocate(part     (   n_verts));  part     (:)   = 1
  allocate(row      ( 1+n_verts))
  allocate(col      ( 2*n_edges))

  allocate(arr_1(n_edges))
  allocate(arr_2(n_edges))

  ! Allocate memory for control parameters for METIS
  n_constrains = 1
  n_parts      = 4   ! how many partitions at each level
  allocate(imbalance   (n_constrains))
  allocate(vert_weights(n_verts * n_constrains))
  allocate(edge_weights(n_edges + n_edges))
  allocate(vert_data   (n_verts))
  allocate(part_weight (n_parts * n_constrains))

  do c = 1, grid % n_cells
    grid % level(0) % cell(c) = c
  end do
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
      n_verts = i

      ! Form edge connectivity in current partition at current level
      edges_v(:,:) = 0
      i = 0
      do s = 1, grid % n_faces
        c1 = grid % faces_c(1, s)
        c2 = grid % faces_c(2, s)

        if(c2 > 0) then
          if(grid % level(sl+1) % cell(c1) .eq. lp .and.  &
             grid % level(sl+1) % cell(c2) .eq. lp) then
            i = i + 1
            edges_v(1,i) = new_c(c1)
            edges_v(2,i) = new_c(c2)
          end if
        end if
      end do
      n_edges = i

      ! Form stars (cell to cell connectivity) for partition and level
      star_size(:) = 0
      star   (:,:) = 0
      do s = 1, n_edges
        nc1 = edges_v(1,s)
        nc2 = edges_v(2,s)

        star_size(nc1) = star_size(nc1) + 1
        star_size(nc2) = star_size(nc2) + 1

        star(star_size(nc1), nc1) = nc2
        star(star_size(nc2), nc2) = nc1
      end do

      ! Fill-up the structures needed to call METIS (row and col)

      row(:) = 0
      row(1) = 0
      do nc = 1, n_verts
        row(nc+1) = row(nc) + star_size(nc)
      end do

      col(:) = 0
      do nc = 1, n_verts
        do i = 1, star_size(nc)
          col(row(nc) + i) = star(i, nc) - 1  ! -1, METIS works from 0
        end do
      end do

      imbalance(:)    = 1.001
      vert_weights(:) = 1
      edge_weights(:) = 1
      vert_data(:)    = 1
      part_weight(:)  = 1.0 / real(n_parts)

      !-----------------------------------------!
      !   Exectute the call to METIS function   !
      !-----------------------------------------!

      metis_options = -1                       ! Initialize all to default
      metis_options(METIS_OPTION_DBGLVL) = 0
      metis_options(METIS_OPTION_CONTIG) = 1

      call METIS_PartGraphRecursive(n_verts,       &  !  1. (in), int
                                    n_constrains,  &  !  2. (in), int
                                    row,           &  !  3. (in), int(:)
                                    col,           &  !  4. (in), int(:)
                                    vert_weights,  &  !  5. (in), int(:)
                                    vert_data,     &  !  6. (in), int(:)
                                    edge_weights,  &  !  7. (in), int(:)
                                    n_parts,       &  !  8. (in), int(:)
                                    part_weight,   &  !  9. (in), real(:)
                                    imbalance,     &  ! 10. (in), real(:)
                                    metis_options, &  ! 11. (in), int(:)
                                    return_val,    &  ! 12. (out) int(:)
                                    part(:))          ! 13. (out) int(:)

      part(1:n_verts) = part(1:n_verts) + 1  ! +1, METIS works from zero

      !-----------------------------------------------------!
      !   Save the result from the call to METIS function   !
      !-----------------------------------------------------!
      if(lev .eq. 1) then
        do c = 1, n_verts
          grid % level(sl) % cell(c) = part(c)
        end do
      else
        do nc = 1, n_verts
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
    arr_1(:) = 0
    arr_2(:) = 0
    do s = 1, grid % level(lev) % n_faces
      arr_1(s) = grid % level(lev) % faces_c(1, s)
      arr_2(s) = grid % level(lev) % faces_c(2, s)
    end do

    ! Sort faces by first cell index
    call Sort_Mod_Int_Carry_Int(arr_1(1:grid % level(lev) % n_faces),  &
                                arr_2(1:grid % level(lev) % n_faces))

    write(lev+100,*), 'Level: ', lev, 'after sorting faces by first index'
    do s = 1, grid % level(lev) % n_faces
      write(lev+100,'(3i6)'), s, arr_1(s), arr_2(s)
    end do

    ! Sort faces by second cell index
    arr_s = 1
    val_1 = arr_1(arr_s)
    do s = 2, grid % level(lev) % n_faces
      if(arr_1(s) .ne. val_1) then
        arr_e = s - 1
        if(arr_e > arr_s) then
          call Sort_Mod_Int_Carry_Int(arr_1(arr_s:arr_e),  &
                                      arr_2(arr_s:arr_e))
        end if
        arr_s = arr_e + 1
        val_1 = arr_1(arr_s)
      end if
    end do

    write(lev+200,*), 'Level: ', lev, 'after sorting faces by second index'
    do s = 1, grid % level(lev) % n_faces
      write(lev+200,'(3i6)'), s, arr_1(s), arr_2(s)
    end do

    ! Selece faces by both cell indices
    i     = 1
    val_1 = arr_1(1)
    val_2 = arr_2(1)
    grid % level(lev) % faces_c(1, i) = val_1
    grid % level(lev) % faces_c(2, i) = val_2

    do s = 2, grid % level(lev) % n_faces
      if(arr_1(s) .ne. val_1 .or.  &
         arr_2(s) .ne. val_2) then
        i = i + 1
        grid % level(lev) % faces_c(1, i) = arr_1(s)
        grid % level(lev) % faces_c(2, i) = arr_2(s)
        val_1 = arr_1(s)
        val_2 = arr_2(s)
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

  deallocate(col)
  deallocate(row)
  deallocate(part_weight)
  deallocate(vert_data)
  deallocate(edge_weights)
  deallocate(vert_weights)
  deallocate(imbalance)
  deallocate(old_c)
  deallocate(new_c)

  end subroutine
