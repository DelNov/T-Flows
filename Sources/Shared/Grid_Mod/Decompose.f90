!==============================================================================!
  subroutine Grid_Mod_Decompose(grid, n_parts)
!------------------------------------------------------------------------------!
!   Coarsens the grid with METIS library.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n_parts
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, s, i
  integer              :: n_verts, n_edges
  integer, allocatable :: star(:,:), star_size(:),  &
                          edges_v(:,:), edges_c(:)

! Variabes for call to METIS
  integer              :: n_constrains,      &  ! number of constraints
                          return_val            ! return value from METIS
  integer, allocatable :: row(:),            &  ! incidency matrix in ...
                          col(:),            &  ! compresser row format
                          vert_weights(:),   &  ! weights of vertices
                          edge_weights(:),   &  ! weights of edges
                          vert_data(:),      &  ! amount of data for vertices
                          part(:)               ! result of partitioning
  real, allocatable    :: part_weight(:),  &
                          imbalance(:)          ! allowed imbalance
  integer              :: metis_options(41)     ! options passed to METIS
!==============================================================================!

  ! Number of vertices and number of edges for first level
  n_verts = grid % n_cells
  n_edges = 0
  do s = 1, grid % n_faces
    c2 = grid % faces_c(2,s)
    if(grid % faces_c(2,s) > 0) n_edges = n_edges + 1
  end do

  ! Once n_verts(1) and n_edegs(1) are known, allocate memory
  allocate(edges_v  ( 2, n_edges));  edges_v  (:,:) = 0
  allocate(edges_c  (    n_edges));  edges_c  (:)   = 0
  allocate(star_size(    n_verts));  star_size(:)   = 0
  allocate(star     (24, n_verts));  star     (:,:) = 0
  allocate(part     (    n_verts));  part     (:)   = 0

  ! Form edge connectivity
  i = 0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then
      i = i + 1
      edges_v(1,i) = c1
      edges_v(2,i) = c2
    end if
  end do

  !------------------------------!
  !   Form stars at this level   !
  !------------------------------!
  do s = 1, n_edges
    c1 = edges_v(1,s)
    c2 = edges_v(2,s)

    star_size(c1) = star_size(c1) + 1
    star_size(c2) = star_size(c2) + 1

    star(star_size(c1), c1) = c2
    star(star_size(c2), c2) = c1
  end do

  !---------------------------------------------------------------!
  !   Fill-up the structures needed to call METIS (row and col)   !
  !---------------------------------------------------------------!

  ! Fill up the rows
  allocate(row(n_verts + 1))
  row(1) = 0
  do c = 1, n_verts
    row(c+1) = row(c) + star_size(c)
  end do

  ! Fill up columns
  allocate(col(row(n_verts + 1)))
  do c = 1, n_verts
    do i = 1, star_size(c)
      col(row(c) + i) = star(i, c) - 1  ! -1, METIS works from 0
    end do
  end do

  n_constrains =  1
  allocate(imbalance (n_constrains));            imbalance(:)    = 1.001
  allocate(vert_weights(n_verts*n_constrains));  vert_weights(:) = 1
  allocate(edge_weights(row(n_verts+1)));        edge_weights(:) = 1
  allocate(vert_data(n_verts));                  vert_data(:)    = 1
  allocate(part_weight(n_parts*n_constrains));   part_weight(:)  &
                                                      = 1.0 / real(n_parts)
  !-----------------------------------------!
  !   Exectute the call to METIS function   !
  !-----------------------------------------!

  metis_options = -1                       ! Initialize all to default
  metis_options(METIS_OPTION_DBGLVL) = 0

  call Metis_PartGraphRecursive(n_verts,       &  !  1. (in), int
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

  part(:) = part(:) + 1  ! +1, METIS works from zero

  deallocate(part_weight)
  deallocate(vert_data)
  deallocate(edge_weights)
  deallocate(vert_weights)
  deallocate(imbalance)
  deallocate(col)
  deallocate(row)

  !-----------------------------------------------------!
  !   Save the result from the call to METIS function   !
  !-----------------------------------------------------!
  do c = 1, grid % n_cells
    grid % comm % cell_proc(c) = part(c)
  end do

  end subroutine
