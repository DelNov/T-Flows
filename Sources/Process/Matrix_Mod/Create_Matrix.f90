!==============================================================================!
  subroutine Create_Matrix(Matrix, Grid)
!------------------------------------------------------------------------------!
!   Determines the topology of the system matrix.                              !
!   It relies only on faces_c structure. Try to keep it that way.              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Matrix_Type)       :: Matrix
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, pos, pos1, pos2, n
  integer              :: c1, c2  ! cell 1 and 2
  integer              :: n1, n2  ! neighbour 1 and 2
  integer, allocatable :: stencw(:)
!==============================================================================!

  ! Store pointer to the grid
  Matrix % pnt_grid => Grid

  ! Allocate memory for matrix
  allocate(Matrix % row(  Grid % n_cells + 1));                Matrix % row=0
  allocate(Matrix % dia(  Grid % n_cells));                    Matrix % dia=0
  allocate(Matrix % fc (  Grid % n_faces));                    Matrix % fc =0.0
  allocate(Matrix % sav( -Grid % n_bnd_cells:Grid % n_cells)); Matrix % sav=0.0
  allocate(Matrix % pos(2,Grid % n_faces));                    Matrix % pos=0

  ! Allocate memory for local array
  allocate(stencw(Grid % n_cells)); stencw = 1

  ! Initialize number of nonzeros
  Matrix % nonzeros = 0

  ! Compute stencis widths
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  ! Count the nonzero entries and allocate the memory for the array 
  n = 0
  do c = 1, Grid % n_cells
    n = n + stencw(c)
  end do
  Matrix % nonzeros = n + 1
  allocate(Matrix % val(n+1)); Matrix % val = 0. ! it reffers to Matrix % row+1 
  allocate(Matrix % col(n+1)); Matrix % col = 0  ! it reffers to Matrix % row+1 
  allocate(Matrix % mir(n+1)); Matrix % mir = 0  ! it reffers to Matrix % row+1 

  ! Form Matrix % row and diagonal only formation of Matrix % col
  Matrix % row(1) = 1
  do c = 1, Grid % n_cells
    Matrix % row(c + 1) = Matrix % row(c) + stencw(c)
    Matrix % col(Matrix % row(c)) = c   ! it is first to its own self
    stencw(c) = 1                       ! reset stencw
  end do 

  ! Extend Matrix % col entries with off-diagonal terms
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      Matrix % col(Matrix % row(c1) + stencw(c1)) = c2
      Matrix % col(Matrix % row(c2) + stencw(c2)) = c1
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  ! Sort Matrix % col to make them nice and neat
  ! and also locate the position of diagonal
  do c = 1, Grid % n_cells
    call Sort % Int_Array(Matrix % col(Matrix % row(c) :  &
                          Matrix % row(c) + stencw(c) - 1))
    do pos = Matrix % row(c), Matrix % row(c+1)-1
      if(Matrix % col(pos) .eq. c) Matrix % dia(c) = pos
    end do
  end do 

  ! Find mirror positions
  do c1 = 1, Grid % n_cells
    do pos1 = Matrix % row(c1), Matrix % row(c1 + 1) - 1
      n1 = Matrix % col(pos1)  ! at this point you have c1 and n1

      ! Inner loop (it might probably go from 1 to c1-1
      c2 = n1
      do pos2 = Matrix % row(c2), Matrix % row(c2 + 1) - 1
        n2 = Matrix % col(pos2)  ! at this point you have c2 and n2

        if(n2 .eq. c1) then
          Matrix % mir(pos1) = pos2
          Matrix % mir(pos2) = pos1
          goto 2  ! done with the inner loop, get out
        end if
      end do

2     continue
    end do
  end do

  ! Connect faces with matrix entries
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then

      ! Where is matrix(c1,c2) and ...
      do c = Matrix % row(c1), Matrix % row(c1+1)-1 
        if(Matrix % col(c) .eq. c2) Matrix % pos(1, s) = c
      end do

      ! ... where is matrix(c2,c1)
      do c=Matrix % row(c2),Matrix % row(c2+1)-1 
        if(Matrix % col(c) .eq. c1) Matrix % pos(2, s) = c
      end do
    end if
  end do

  ! Bare-bone coefficients for system matrix with overrelaxed correction
  ! In Denner's thesis at page 51 it is used: alpha_f = 1.0 / (n dot d)
  ! (d is s in Denner's thesis) or: alpha_f = sqrt(s dot s) / (s dot d).
  ! If multiplied with surface area which is sqrt(s dot s), the final
  ! coefficient reads: fc = (s dot s) / (s dot d), as computed below:
  do s = 1, Grid % n_faces
    Matrix % fc(s) = (  Grid % sx(s) * Grid % sx(s)    &
                      + Grid % sy(s) * Grid % sy(s)    &
                      + Grid % sz(s) * Grid % sz(s) )  &
                   / (  Grid % dx(s) * Grid % sx(s)    &
                      + Grid % dy(s) * Grid % sy(s)    &
                      + Grid % dz(s) * Grid % sz(s) )
  end do

  deallocate(stencw)

  end subroutine
