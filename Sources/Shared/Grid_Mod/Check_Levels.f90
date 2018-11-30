!==============================================================================!
  subroutine Grid_Mod_Check_Levels(grid)
!------------------------------------------------------------------------------!
!   Coarsens the grid with METIS library.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, c1, c2, nc, nc1, nc2, s, i, lev, lev_parts
  integer              :: arr_s, arr_e, val_1, val_2
  integer              :: c_lev, c1_lev_c, c2_lev_c, s_lev
  integer, allocatable :: cell_mapping(:,:)
  integer, allocatable :: face_hits(:)
!==============================================================================!

  print *, '#====================================='
  print *, '# Checking sanity of multigrid levels '
  print *, '#-------------------------------------'

  !--------------------!
  !                    !
  !   Sanity check 1   !
  !                    !
  !--------------------!
  print *, '# Sanity check 1 '
  i = 0
  do lev = 1, grid % n_levels
    i = max(i, grid % level(lev) % n_cells)
  end do
  allocate(cell_mapping(MAX_MG_LEVELS, i)); cell_mapping = 0

  do lev = 1, grid % n_levels - 1
    print '(a,i2,a,i2)', ' # ... levels', lev, ' and', lev+1

    ! Browse through parts of this level
    do c_lev = 1, grid % level(lev) % n_cells
      do c = 1, grid % n_cells
        if(grid % level(lev) % cell(c) == c_lev) then
          if(cell_mapping(lev, c_lev) .eq. 0) then
            cell_mapping(lev, c_lev) = grid % level(lev+1) % cell(c)
          else
            if(cell_mapping(lev, c_lev) .ne. grid % level(lev+1) % cell(c)) then
              print *, '# Mapping failed at level ', lev
              print *, '# Stopping the program!   '
              stop
            end if
          end if
        end if
      end do
    end do

  end do  ! lev

  !--------------------!
  !                    !
  !   Sanity check 2   !
  !                    !
  !--------------------!
  print *, '# Sanity check 2 '
  do lev = 1, grid % n_levels
    do s = 1, grid % level(lev) % n_faces
      c1 = grid % level(1) % faces_c(1, s)
      c2 = grid % level(1) % faces_c(2, s)

      c1_lev_c = min(grid % level(lev) % cell(c1),  &
                     grid % level(lev) % cell(c2))
      c2_lev_c = max(grid % level(lev) % cell(c1),  &
                     grid % level(lev) % cell(c2))

      s_lev = grid % level(lev) % face(s)

      if(s_lev > 0) then
        if(grid % level(lev) % faces_c(1,s_lev) .ne. c1_lev_c .or. &
           grid % level(lev) % faces_c(2,s_lev) .ne. c2_lev_c) then
          print *, '# Cell-face mapping failed at level ', lev
          print *, '# Stopping the program!   '
          stop
        end if
      end if

    end do
  end do

  !--------------------!
  !                    !
  !   Sanity check 3   !
  !                    !
  !--------------------!
  print *, '# Sanity check 3 '
  do lev = 1, grid % n_levels
    do c_lev = 1, grid % level(lev) % n_cells
      do s_lev = 1, grid % level(lev) % n_faces
        if(grid % level(lev) % faces_c(1, s_lev) .eq. c_lev .or.  &
           grid % level(lev) % faces_c(2, s_lev) .eq. c_lev) then
          goto 1
        end if
      end do
      print *, '# Cell ', c_lev, ' at level ', lev, ' not found among the faces'
      print *, '# Stopping the program!   '
      stop
1     continue
    end do

  end do

  !--------------------!
  !                    !
  !   Sanity check 4   !
  !                    !
  !--------------------!
  print *, '# Sanity check 4 '
  do lev = 1, grid % n_levels
    allocate(face_hits(grid % level(lev) % n_faces))
    face_hits(:) = 0
    do s = 1, grid % level(1) % n_faces
      s_lev = grid % level(lev) % face(s)
      if(s_lev > 0) then
        face_hits(s_lev) = face_hits(s_lev) + 1
      end if
    end do
    do s_lev = 1, grid % level(lev) % n_faces
      if(face_hits(s_lev) .eq. 0) then
        print *, '# Sad!  Face ', s_lev, ' on level ', lev, ' was never hit'
        print *, '# Stopping the program!   '
        stop
      end if
    end do
    deallocate(face_hits)
  end do

  end subroutine
