!==============================================================================!
  subroutine Create_Sparse(A, Grid)
!------------------------------------------------------------------------------!
!   Note:                                                                      !
!                                                                              !
!   * This subroutine is executed on CPU (host) only.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Type)      :: A   !! parent connectivity matrix class
  type(Grid_Type), target :: Grid  !! grid on which it is created
!-----------------------------------[Locals]-----------------------------------!
  integer :: non_z, non_z_glo, run
  integer :: c, i_cel, s, c1, c2
  integer :: col_a, row_a, pos_a
!==============================================================================!

  ! Store pointer to the grid
  A % pnt_grid => Grid

  O_Print '(a)', ' # Creating a sparse connectivity matrix'

  !--------------------------------!
  !                                !
  !   Form the compressed matrix   !
  !                                !
  !--------------------------------!
  do run = 1, 2  ! first run is just for counting non-zeros, second for filling

    non_z = 0    ! reset the counter of non-zero entries

    !---------------------------------------------------!
    !   Browse through cells and also through buffers   !
    !- - - - - - - - - - - - - - - - - - - - - - - - - -!
    !   (When I was avoiding cells in buffers, I had    !
    !   "holes" in CSR format and a lot of troubles)    !
    !---------------------------------------------------!
    do c1 = Cells_In_Domain_And_Buffers()  ! this whole routine is on CPU

      ! Set this row index, one cell will be present for sure, own self
      if(run .eq. 2) A % row(c1) = non_z + 1

      !-------------------------------------------------------!
      !   Store the central entry, be it in obstacle or not   !
      !-------------------------------------------------------!
      non_z = non_z + 1
      if(run .eq. 2) A % col(non_z) = c1

      do i_cel = 1, Grid % cells_n_cells(c1)
        c2 = Grid % cells_c(i_cel, c1)

        if(c2 .gt. 0) then
          non_z = non_z + 1
          if(run .eq. 2) A % col(non_z) = c2
        end if    ! c2 is inside cell
      end do      ! i_cel

    end do        ! c1

    ! Work out the and print number of non-zeros over all processors, globally
    if(run .eq. 1) then
      non_z_glo = non_z
      call Global % Sum_Int(non_z_glo)
      O_Print '(a,i15)', ' # Number of nonzeros: ', non_z_glo
    end if

    !--------------------------------------------------------------!
    !   If this is the end of the first run, allocate the memory   !
    !--------------------------------------------------------------!
    if(run .eq. 1) then
      A % nonzeros = non_z
      allocate(A % row(Grid % n_cells+1));   A % row = 0
      allocate(A % dia(Grid % n_cells));     A % dia = 0
      allocate(A % col(non_z));              A % col = 0
      allocate(A % fc (Grid % n_faces));     A % fc  = 0
      Assert(Grid % n_faces .gt. 0)
      allocate(A % pos(2, Grid % n_faces));  A % pos = 0
    end if

    !-----------------------------------------------------!
    !   Wrap it up - set the end of the last cell's row   !
    !-----------------------------------------------------!
    if(run .eq. 2) A % row(Grid % n_cells + 1) = non_z + 1

  end do  ! run

  !--------------------------------------!
  !                                      !
  !   Sort each row in ascending order   !
  !                                      !
  !--------------------------------------!
  do c = Cells_In_Domain_And_Buffers()  ! this whole routine is on CPU
    call Sort % Int_Array(A % col(A % row(c) : A % row(c+1)-1))
  end do

  !---------------------------------!
  !                                 !
  !   Find positions of diagonals   !
  !                                 !
  !---------------------------------!
  do row_a = Cells_In_Domain_And_Buffers()  ! this whole routine is on CPU
    do pos_a = A % row(row_a), A % row(row_a + 1) - 1
      col_a = A % col(pos_a)  ! at this point you have row_a and col_a
      if(col_a == row_a) then
        A % dia(row_a) = pos_a
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

    A % fc(s) = (   Grid % sx(s) * Grid % sx(s)    &
                    + Grid % sy(s) * Grid % sy(s)    &
                    + Grid % sz(s) * Grid % sz(s) )  &
                 / (  Grid % dx(s) * Grid % sx(s)    &
                    + Grid % dy(s) * Grid % sy(s)    &
                    + Grid % dz(s) * Grid % sz(s) )
  end do

  !---------------------------------------!
  !                                       !
  !   Connect faces with matrix entries   !
  !                                       !
  !---------------------------------------!
  do s = Faces_In_Domain_And_At_Buffers()  ! this whole routine is on CPU
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Where is matrix(c1,c2) and ...
    do c = A % row(c1), A % row(c1+1)-1
      if(A % col(c) .eq. c2) then
        A % pos(1, s) = c
        exit
      end if
    end do

    ! ... where is matrix(c2,c1)
    do c = A % row(c2),A % row(c2+1)-1
      if(A % col(c) .eq. c1) then
        A % pos(2, s) = c
        exit
      end if
    end do

    ! These connections shouldn't be left at zeroes
    Assert(A % pos(1,s) .ne. 0)
    Assert(A % pos(2,s) .ne. 0)

  end do

  !------------------------------------------------!
  !                                                !
  !   Perform a test: face-based matrix pointers   !
  !   must be inside the bounds defined by "row"   !
  !                                                !
  !------------------------------------------------!
  do s = Faces_In_Domain_And_At_Buffers()  ! this whole routine is on CPU
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    Assert(A % pos(1,s) .ge. A % row(c1))
    Assert(A % pos(1,s) .le. A % row(c1+1)-1)

    if(c2 .gt. 0) then
      Assert(A % pos(2,s) .ge. A % row(c2))
      Assert(A % pos(2,s) .le. A % row(c2+1)-1)
    end if
  end do

  !---------------------------------------------------------!
  !                                                         !
  !    Finally, allocate memory for storing matrix values   !
  !                                                         !
  !---------------------------------------------------------!
  allocate(A % val  (A % nonzeros));     A % val   = 0.0
  allocate(A % d_inv(Grid % n_cells));   A % d_inv = 0.0

  end subroutine
