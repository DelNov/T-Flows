!==============================================================================!
  subroutine Create_Sparse_Con(Con, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Con_Type)   :: Con   !! parent connectivity matrix class
  type(Grid_Type), target :: Grid  !! grid on which it is created
!-----------------------------------[Locals]-----------------------------------!
  integer :: non_z, run
  integer :: c, i_cel, s, c1, c2
  integer :: col_a, row_a, pos_a
!==============================================================================!

  ! Store pointer to the grid
  Con % pnt_grid => Grid

  print '(a)', ' # Creating a sparse matrix'

  !--------------------------------!
  !                                !
  !   Form the compressed matrix   !
  !                                !
  !--------------------------------!
  do run = 1, 2  ! first run is just for counting non-zeros, second for filling

    non_z = 0    ! reset the counter of non-zero entries

    !----------------------------------!
    !   Browse through all the cells   !
    !----------------------------------!
    do c1 = 1, Grid % n_cells

      ! Set this row index, one cell will be present for sure, own self
      if(run .eq. 2) Con % row(c1) = non_z + 1

      !-------------------------------------------------------!
      !   Store the central entry, be it in obstacle or not   !
      !-------------------------------------------------------!
      non_z = non_z + 1
      if(run .eq. 2) Con % col(non_z) = c1

      do i_cel = 1, Grid % cells_n_cells(c1)
        c2 = Grid % cells_c(i_cel, c1)

        if(c2 .gt. 0) then
          non_z = non_z + 1
          if(run .eq. 2) Con % col(non_z) = c2
        end if    ! c2 is inside cell
      end do      ! i_cel

    end do        ! c1

    print '(a,i15)', ' # Number of nonzeros: ', non_z

    !--------------------------------------------------------------!
    !   If this is the end of the first run, allocate the memory   !
    !--------------------------------------------------------------!
    if(run .eq. 1) then
      Con % nonzeros = non_z
      allocate(Con % row(Grid % n_cells+1));   Con % row = 0
      allocate(Con % dia(Grid % n_cells));     Con % dia = 0
      allocate(Con % col(non_z));              Con % col = 0
      allocate(Con % fc (Grid % n_faces));     Con % fc  = 0
      Assert(Grid % n_faces .gt. 0)
      allocate(Con % pos(2, Grid % n_faces));  Con % pos = 0
    end if

    !-----------------------------------------------------!
    !   Wrap it up - set the end of the last cell's row   !
    !-----------------------------------------------------!
    if(run .eq. 2) Con % row(Grid % n_cells + 1) = non_z + 1

  end do  ! run

  !--------------------------------------!
  !                                      !
  !   Sort each row in ascending order   !
  !                                      !
  !--------------------------------------!
  do c = 1, Grid % n_cells
    call Sort % Int_Array(Con % col(Con % row(c) : Con % row(c+1)-1))
  end do

  !---------------------------------!
  !                                 !
  !   Find positions of diagonals   !
  !                                 !
  !---------------------------------!
  do row_a = 1, Grid % n_cells
    do pos_a = Con % row(row_a), Con % row(row_a + 1) - 1
      col_a = Con % col(pos_a)  ! at this point you have row_a and col_a
      if(col_a == row_a) then
        Con % dia(row_a) = pos_a
        goto 1
      end if
    end do
1   continue
  end do

  !-----------------------------------!
  !                                   !
  !   Bare-bone matrix coefficients   !
  !                                   !
  !-----------------------------------!
  do s = 1, Grid % n_faces

    Assert(Grid % s(s) .gt. TINY)
    Assert(Grid % d(s) .gt. TINY)

    Con % fc(s) = Grid % s(s) / Grid % d(s)

  end do

  !---------------------------------------!
  !                                       !
  !   Connect faces with matrix entries   !
  !                                       !
  !---------------------------------------!
  do s = Grid % n_bnd_cells + 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Where is matrix(c1,c2) and ...
    do c = Con % row(c1), Con % row(c1+1)-1
      if(Con % col(c) .eq. c2) then
        Con % pos(1, s) = c
        exit
      end if
    end do

    ! ... where is matrix(c2,c1)
    do c=Con % row(c2),Con % row(c2+1)-1
      if(Con % col(c) .eq. c1) then
        Con % pos(2, s) = c
        exit
      end if
    end do

  end do

  end subroutine
