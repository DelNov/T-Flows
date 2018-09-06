!==============================================================================!
  subroutine Matrix_Mod_Topology(grid, matrix)
!------------------------------------------------------------------------------!
!   Determines the topology of the system matrix.                              !
!   It relies only on SideC structure. Try to keep it that way.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  type(Matrix_Type)       :: matrix
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, pos, pos1, pos2, n
  integer              :: c1, c2  ! cell 1 and 2
  integer              :: n1, n2  ! neighbour 1 and 2
  integer, allocatable :: stencw(:)
!==============================================================================!

  ! Store pointer to the grid
  matrix % pnt_grid => grid

  ! Memory allocation
  allocate(stencw(grid % n_cells)); stencw = 1

  if(this_proc < 2) print *, '# Determining matrix topology.'

  ! Initialize number of nonzeros
  matrix % nonzeros = 0

  ! Compute stencis widths
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  ! Count the nonzero entries and allocate the memory for the array 
  n = 0
  do c = 1, grid % n_cells
    n = n + stencw(c)
  end do
  matrix % nonzeros = n + 1
  allocate(matrix % val(n+1)); matrix % val = 0. ! it reffers to matrix % row+1 
  allocate(matrix % col(n+1)); matrix % col = 0  ! it reffers to matrix % row+1 
  allocate(matrix % mir(n+1)); matrix % mir = 0  ! it reffers to matrix % row+1 

  ! Form matrix % row and diagonal only formation of matrix % col
  matrix % row(1) = 1
  do c = 1, grid % n_cells
    matrix % row(c + 1) = matrix % row(c) + stencw(c)
    matrix % col(matrix % row(c)) = c   ! it is first to its own self
    stencw(c) = 1                       ! reset stencw
  end do 

  ! Extend matrix % col entries with off-diagonal terms
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
      matrix % col(matrix % row(c1) + stencw(c1)) = c2
      matrix % col(matrix % row(c2) + stencw(c2)) = c1
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  ! Sort matrix % col to make them nice and neat
  ! and also locate the position of diagonal
  do c = 1, grid % n_cells
    call Sort_Mod_Int(matrix % col(matrix % row(c) :  &
                                   matrix % row(c) + stencw(c) - 1))
    do pos = matrix % row(c),matrix % row(c+1)-1
      if(matrix % col(pos) .eq. c) matrix % dia(c) = pos
    end do
  end do 

  ! Find mirror positions
  do c1 = 1, grid % n_cells
    do pos1 = matrix % row(c1), matrix % row(c1 + 1) - 1 
      n1 = matrix % col(pos1)  ! at this point you have c1 and n1  

      ! Inner loop (it might probably go from 1 to c1-1
      c2 = n1
      do pos2 = matrix % row(c2), matrix % row(c2 + 1) - 1 
        n2 = matrix % col(pos2)  ! at this point you have c2 and n2 

        if(n2 .eq. c1) then
          matrix % mir(pos1) = pos2 
          matrix % mir(pos2) = pos1 
          goto 2  ! done with the inner loop, get out
        end if
      end do

2     continue
    end do
  end do

  ! Connect faces with matrix entries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then

      ! Where is matrix(c1,c2) and ...
      do c = matrix % row(c1), matrix % row(c1+1)-1 
        if(matrix % col(c)  .eq.  c2) matrix % pos(1, s) = c
      end do

      ! Where is matrix(c2,c1) and ...
      do c=matrix % row(c2),matrix % row(c2+1)-1 
        if(matrix % col(c)  .eq.  c1) matrix % pos(2, s) = c
      end do
    end if
  end do

  if(this_proc < 2) print *, '# Finished !'
 
  deallocate(stencw)

  end subroutine
