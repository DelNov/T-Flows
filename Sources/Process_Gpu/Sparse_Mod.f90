#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Sparse_Mod
!------------------------------------------------------------------------------!
!>  The Sparse_Mod module is designed to handle the storage and manipulation of
!>  sparse matrices within T-Flows. It stores matrices in the compressed row
!>  storage (CSR) format, a common format for storing sparse matrices.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Sparse type                                                                !
!                                                                              !
!   Description:                                                               !
!   - The Sparse_Con_Type and Sparse_Val_Type are  custom data type for        !
!     handling sparse matrices in T-Flows.  They uses the compressed row       !
!     storage (CRS) format, an efficient method for storing sparse matrices.   !
!     (See: http://netlib.org/linalg/html_templates/node91.html)               !
!                                                                              !
!   - The Sparse_Con_Type holds only the fields describing the connectivity    !
!     in a matrix, and Sparse_Val_Type stores the values.  Clearly, the        !
!     latter is variable dependent (it will not be the same for momentum and   !
!     pressure, for example, but the former will.  Storing them separatelly    !
!     leads to memory savings, particularly important on GPUs.                 !
!                                                                              !
!   Compressed Row Format:                                                     !
!   - In CRS, a matrix is represented by three arrays: 'val' for non-zero      !
!     values, 'col' for column indices, and 'row' for row pointers. This       !
!     format is advantageous for matrix-vector multiplication and other        !
!     matrix operations.                                                       !
!                                                                              !
!   Example:                                                                   !
!   - Illustrates a 4x4 matrix and its CRS representation.  The 'dia' array    !
!     stores diagonal positions, and 'pos' holds face-based positions in the   !
!     matrix. 'glo' represents global cell numbering, differing from T-Flows'  !
!     internal numbering, but compatible with PETSc.                           !
!                                                                              !
!       c   c  .    c                                                          !
!       o   o  .    o                                                          !
!       l   l       l                                                          !
!                                                                              !
!       1   2       n                                                          !
!                                                                              !
!     [ 10   0   4   5 ]  --> row 1                                            !
!     [  2  12  -1   0 ]  --> rows store discretized control volumes           !
!     [  0   1  99   7 ]  ...                                                  !
!     [ -3  11   0  53 ]  --> row n                                            !
!                                                                              !
!     Compressed row storage of the above matrix reads:                        !
!                                                                              !
!     A % val = [  10   4   5   2  12  -1   1  99   7  -3  11  53 ]            !
!     A % col = [   1   3   4   1   2   3   2   3   4   1   2   4 ]            !
!     A % row = [   1   4   7  10 ]                                            !
!                                                                              !
!     A % dia = [   1   5   9  12 ]                                            !
!                                                                              !
!==============================================================================!

  !------------------------------!
  !                              !
  !   Sparse connectivity type   !
  !                              !
  !------------------------------!
  type Sparse_Con_Type
  !> A type which encapsulates matrix connectivity in CRS format.

    type(Grid_Type), pointer :: pnt_grid  !! pointer to grid

    integer              :: nonzeros  !! number of nonzero entries
    real,    allocatable :: fc (:)    !! bare matrix entry for face
    integer, allocatable :: col(:)    !! column positions
    integer, allocatable :: row(:)    !! beginning of each row
    integer, allocatable :: dia(:)    !! diagonal positions
    integer, allocatable :: pos(:,:)  !! face-based position of the matrix

    contains
      procedure :: Create_Sparse_Con

  end type

  !------------------------!
  !                        !
  !   Sparse values type   !
  !                        !
  !------------------------!
  type Sparse_Val_Type
  !> A type which encapsulates matrix values in CRS format.

    type(Grid_Type), pointer :: pnt_grid  !! pointer to grid

    logical :: formed = .false.  !! set to true when matrix is formed

    real, allocatable :: val(:)    !! value
    real, allocatable :: d_inv(:)  !! inverse diagonal for preconditioner

    contains
      procedure :: Create_Sparse_Val

  end type

  contains
#   include "Sparse_Mod/Create_Sparse_Con.f90"
#   include "Sparse_Mod/Create_Sparse_Val.f90"

  end module
