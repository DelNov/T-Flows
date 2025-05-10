!==============================================================================!
  subroutine Create_Matrix(A, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine serves as a constructor of Matrix_Type and it allocates
!>  memory for storage of system matrices in CRS formant and is establishes
!>  the connectivity of the matrix based on the computational grid.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Grid linkage: Establishes a connection to the computational grid.        !
!   * Memory allocation: Allocates memory for matrix components such as        !
!     row indices, column indices, and non-zero values.                        !
!   * Stencil width calculation: Determines the width of each stencil based on !
!     the grid's face-cell connectivity.                                       !
!   * Assembly of CRS format: Constructs the compressed row storage format for !
!     the matrix, including off-diagonal terms.                                !
!   * Sorting: Orders the column indices to maintain the structure of the CRS  !
!     format and identifies diagonal positions.                                !
!   * Face-matrix connection: Links faces to corresponding matrix entries.     !
!   * Bare-bone coefficient calculation: Computes preliminary coefficients     !
!     for the system matrix.                                                   !
!   * PETSc global numbering: Assigns global numbering for PETSc, which starts !
!     from zero and differs from T-Flows' internal numbering.                  !
!   * Parallel considerations: Handles global numbering in parallel runs,      !
!     ensuring consistency across different processors.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Matrix_Type)      :: A     !! matrix to be created
  type(Grid_Type), target :: Grid  !! computational grid for matrix assembly
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, pos, n
  integer              :: c1, c2  ! cell 1 and 2
  integer              :: m_lower, m_upper, start
  integer, allocatable :: stencw(:)
  integer, allocatable :: all_lower_ms(:)
!==============================================================================!

  ! Store pointer to the grid
  A % pnt_grid => Grid

  !--------------------------------!
  !   Allocate memory for matrix   !
  !--------------------------------!
  allocate(A % row(  Grid % n_cells + 1));                A % row=0
  allocate(A % dia(  Grid % n_cells));                    A % dia=0
  allocate(A % fc (  Grid % n_faces));                    A % fc =0.0
  allocate(A % v_m( -Grid % n_bnd_cells:Grid % n_cells)); A % v_m=0.0
  allocate(A % pos(2,Grid % n_faces));                    A % pos=0

  !-------------------------------------!
  !   Allocate memory for local array   !
  !-------------------------------------!
  allocate(stencw(Grid % n_cells)); stencw = 1

  !-----------------------------------!
  !   Initialize number of nonzeros   !
  !-----------------------------------!
  A % nonzeros = 0

  !----------------------------!
  !   Compute stencis widths   !
  !----------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  !---------------------------------------------------------------------!
  !   Count the nonzero entries and allocate the memory for the array   !
  !---------------------------------------------------------------------!
  n = 0
  do c = Cells_In_Domain_And_Buffers()
    n = n + stencw(c)
  end do
  A % nonzeros = n + 1
  allocate(A % val(n+1)); A % val = 0. ! it reffers to A % row+1
  allocate(A % col(n+1)); A % col = 0  ! it reffers to A % row+1

  !---------------------------------------------------------!
  !   Form A % row and diagonal only formation of A % col   !
  !---------------------------------------------------------!
  A % row(1) = 1
  do c = Cells_In_Domain_And_Buffers()
    A % row(c + 1) = A % row(c) + stencw(c)
    A % col(A % row(c)) = c   ! it is first to its own self
    stencw(c) = 1                       ! reset stencw
  end do

  !----------------------------------------------------!
  !   Extend A % col entries with off-diagonal terms   !
  !----------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      A % col(A % row(c1) + stencw(c1)) = c2
      A % col(A % row(c2) + stencw(c2)) = c1
      stencw(c1) = stencw(c1) + 1
      stencw(c2) = stencw(c2) + 1
    end if
  end do

  !----------------------------------------------!
  !   Sort A % col to make them nice and neat    !
  !   and also locate the position of diagonal   !
  !----------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
    call Sort % Int_Array(A % col(A % row(c) :  &
                          A % row(c) + stencw(c) - 1))
    do pos = A % row(c), A % row(c+1)-1
      if(A % col(pos) .eq. c) A % dia(c) = pos
    end do
  end do

  !---------------------------------------!
  !   Connect faces with matrix entries   !
  !---------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then

      ! Where is matrix(c1,c2) and ...
      do c = A % row(c1), A % row(c1+1)-1 
        if(A % col(c) .eq. c2) A % pos(1, s) = c
      end do

      ! ... where is matrix(c2,c1)
      do c=A % row(c2),A % row(c2+1)-1 
        if(A % col(c) .eq. c1) A % pos(2, s) = c
      end do
    end if
  end do

  !--------------------------------------------------------------------------!
  !   Bare-bone coefficients for system matrix with overrelaxed correction   !
  !   In Denner's thesis at page 51 it is used: alpha_f = 1.0 / (n dot d)    !
  !   (d is s in Denner's thesis) or: alpha_f = sqrt(s dot s) / (s dot d).   !
  !   If multiplied with surface area which is sqrt(s dot s), the final      !
  !   coefficient reads: fc = (s dot s) / (s dot d), as computed below:      !
  !--------------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    A % fc(s) = (   Grid % sx(s) * Grid % sx(s)    &
                  + Grid % sy(s) * Grid % sy(s)    &
                  + Grid % sz(s) * Grid % sz(s) )  &
               / (  Grid % dx(s) * Grid % sx(s)    &
                  + Grid % dy(s) * Grid % sy(s)    &
                  + Grid % dz(s) * Grid % sz(s) )
  end do

  !----------------------------------------+
  !    Create global numbering for PETSc   !
  !----------------------------------------+------------------!
  !    This one has little to do with global numbering from   !
  !    T-Flows and is unique for each number of processors    !
  !-----------------------------------------------------------!

  ! Total number of unknowns and unknowns in this processor only
  m_upper = Grid % Comm % nc_tot
  m_lower = Grid % n_cells - Grid % Comm % n_buff_cells

  ! Dimensions must spread from all boundary cells through all ...
  ! ... buffers cells to successfully use Grid % Exchange_Cells_Int
  allocate(A % glo(-Grid % n_bnd_cells:Grid % n_cells))
  A % glo(:) = 0

  if(Sequential_Run()) then
    A % glo(1:Grid % n_cells)  &
      = Grid % Comm % cell_glo(1:Grid % n_cells) - 1
  else
    start = 1  ! first row
    allocate(all_lower_ms(N_Procs()));  ! allocate array for all m_lowers
    all_lower_ms(:) = 0                 ! important to initialize to zero

    ! Distribute m_lowers among all processors
    all_lower_ms(This_Proc()) = m_lower
    call Global % Sum_Int_Array(N_Procs(), all_lower_ms)

    start = sum(all_lower_ms(1:This_Proc())) - m_lower

    ! Distribute global numbers over other processors
    do c = 1, m_lower
      A % glo(c) = c + start - 1
    end do
    call Grid % Exchange_Cells_Int(A % glo)

  end if  ! sequential or parallel run?

  end subroutine
