!==============================================================================!
  subroutine Find_Faces(grid)
!------------------------------------------------------------------------------!
!   Find faces inside the domain.  To be more specific, it determines:         !
!                                                                              !
!   grid % n_faces       - final number of faces (boundary + inside)           !
!   grid % faces_n_nodes - number of nodes for each face                       !
!   grid % faces_n       - nodes of each face                                  !
!   grid % faces_c       - pair of cells surrounding each face                 !
!                                                                              !
!   Note: boundary faces have been determined in "Grid_Topology"               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer             :: c, c1, c2, n1, n2, n3, f_nod(4), n_f_nod
  integer             :: n_match, j, match_nodes(-1:8) 
  integer             :: i1, i2, k, cnt
  integer             :: fn(6,4)
  real   ,allocatable :: face_coor(:)
  integer,allocatable :: face_cell(:), starts(:), ends(:)
  real                :: very_big
!==============================================================================!

  very_big = max(grid % n_nodes,grid % n_cells)

  allocate(face_coor(grid % n_cells*6));  face_coor = grid % n_nodes * HUGE
  allocate(face_cell(grid % n_cells*6));  face_cell = 0    
  allocate(starts   (grid % n_cells*6));  starts    = 0    
  allocate(ends     (grid % n_cells*6));  ends      = 0    

  !---------------------------------------------------!
  !   Fill the generic coordinates with some values   !
  !---------------------------------------------------!
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) .eq. 4) fn = neu_tet
    if(grid % cells_n_nodes(c) .eq. 5) fn = neu_pyr
    if(grid % cells_n_nodes(c) .eq. 6) fn = neu_wed
    if(grid % cells_n_nodes(c) .eq. 8) fn = neu_hex 
    do j = 1, 6
      if(grid % cells_bnd_color(j,c) .eq. 0) then

        n_f_nod = 0
        f_nod = -1
        do k = 1, 4
          if(fn(j,k) > 0) then
            f_nod(k) = grid % cells_n(fn(j,k), c)
            n_f_nod = n_f_nod + 1
          end if
        end do

        if( n_f_nod >  0 ) then
          if(f_nod(4) > 0) then
            face_coor((c-1)*6+j) =   &
               very_big*(max(f_nod(1), f_nod(2), f_nod(3), f_nod(4)))   &
            +            min(f_nod(1), f_nod(2), f_nod(3), f_nod(4))
          else
            face_coor((c-1)*6+j) =   &
              very_big*(max(f_nod(1), f_nod(2), f_nod(3)))   &
           +            min(f_nod(1), f_nod(2), f_nod(3))
           end if
          face_cell((c-1)*6+j) = c 
        end if 
      end if
    end do
  end do

  !--------------------------------------------------!
  !   Sort the cell faces according to coordinares   !
  !--------------------------------------------------!
  call Sort_Mod_Real_Carry_Int(face_coor(1:grid % n_cells*6),  & 
                               face_cell(1:grid % n_cells*6))

  !------------------------------------------------!
  !   Anotate cell faces with same coordinates     !
  !   (I am afraid that this might be influenced   !
  !      by the numerical round-off errors)        !
  !------------------------------------------------!
  cnt = 1
  starts(1) = 1
  do c=2,grid % n_cells*6
    if( face_coor(c) .ne. face_coor(c-1) ) then
      cnt = cnt + 1
      starts(cnt) = c
      ends(cnt-1) = c-1
    end if
  end do

  !-------------------------------------------!
  !                                           !
  !   Main loop to fill the SideC structure   !
  !                                           !
  !-------------------------------------------!
  do n3 = 1, cnt
    if(starts(n3) .ne. ends(n3)) then
      do i1=starts(n3),ends(n3)
        do i2=i1+1,ends(n3)
          c1 = min(face_cell(i1),face_cell(i2))
          c2 = max(face_cell(i1),face_cell(i2))
          if(c1 .ne. c2) then

            !------------------------------!
            !   Number of matching nodes   !
            !------------------------------!
            n_match     = 0
            match_nodes = 0 
            do n1 = 1, grid % cells_n_nodes(c1)
              do n2 = 1, grid % cells_n_nodes(c2)
                if(grid % cells_n(n1,c1) .eq. grid % cells_n(n2,c2)) then
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
              if(grid % cells_n_nodes(c1) .eq. 4) fn = neu_tet
              if(grid % cells_n_nodes(c1) .eq. 5) fn = neu_pyr
              if(grid % cells_n_nodes(c1) .eq. 6) fn = neu_wed
              if(grid % cells_n_nodes(c1) .eq. 8) fn = neu_hex
              do j = 1, 6
                if(grid % cells_c(j, c1) .eq. 0  .and.   & ! not set yet
                    ( max( match_nodes(fn(j,1)),0 ) + &
                      max( match_nodes(fn(j,2)),0 ) + &
                      max( match_nodes(fn(j,3)),0 ) + &
                      max( match_nodes(fn(j,4)),0 ) .eq. n_match ) ) then
                  grid % n_faces = grid % n_faces + 1 
                  grid % faces_c(1,grid % n_faces) = c1
                  grid % faces_c(2,grid % n_faces) = c2
                  grid % faces_n_nodes(grid % n_faces) = n_match 
                  do k = 1, 4
                    if(fn(j,k) > 0) then
                      grid % faces_n(k,grid % n_faces) =  &
                      grid % cells_n(fn(j,k), c1)                
                    end if
                  end do
                  grid % cells_c(j, c1) = 1 !  -> means: set
                end if
              end do
            end if   ! n_match .ne. 2
          end if   ! c1 .ne. c2
        end do   ! i2
      end do   ! i1
    end if
  end do    ! do n3

  print '(a38,i9)', '# Number of faces:                   ', grid % n_faces

  end subroutine
