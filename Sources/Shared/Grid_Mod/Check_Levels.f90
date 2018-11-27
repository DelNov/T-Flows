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
  integer              :: n_cells, arr_s, arr_e, val_1, val_2
  integer              :: c_lev, c1_lev, c2_lev, s_lev
  integer, allocatable :: cell_mapping(:,:)
!==============================================================================!

  print *, '#====================================='
  print *, '# Checking sanity of multigrid levels '
  print *, '#-------------------------------------'

  !--------------------!
  !                    !
  !   Sanity check 1   !
  !                    !
  !--------------------!
  i = 0
  do lev = 1, grid % n_levels
    i = max(i, grid % level(lev) % n_cells)
  end do
  allocate(cell_mapping(MAX_MG_LEVELS, i)); cell_mapping = 0

  do lev = 1, grid % n_levels - 1
    n_cells = maxval(grid % level(lev) % cell(:))
    print '(a,i2,a,i2)', ' # Checking levels', lev, ' and', lev+1

    ! Browse through parts of this level
    do c_lev = 1, n_cells
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
  do lev = 1, grid % n_levels
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1, s)
      c2 = grid % faces_c(2, s)

      if(c2 > 0) then
        c1_lev = grid % level(lev) % cell(c1)
        c2_lev = grid % level(lev) % cell(c2)
        s_lev  = grid % level(lev) % face(s)

        if(s_lev > 0) then
          if(grid % level(lev) % faces_c(1,s_lev) .ne.  &
             min(c1_lev, c2_lev) .or.                   &
             grid % level(lev) % faces_c(2,s_lev) .ne.  &
             max(c1_lev, c2_lev)) then
            print *, '# Cell-face mapping failed at level ', lev
            print *, '# Stopping the program!   '
            stop
          end if
        end if
      end if

    end do
  end do

  end subroutine
