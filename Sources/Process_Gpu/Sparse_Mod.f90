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
  use Sort_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !--------------------------------------------------------------------!
  !   Sparse type                                                      !
  !                                                                    !
  !   Description:                                                     !
  !   - The Matrix_Type is a custom data type for handling sparse      !
  !     matrices in T-Flows. It uses the compressed row storage (CRS)  !
  !     format, an efficient method for storing sparse matrices.       !
  !     (See: http://netlib.org/linalg/html_templates/node91.html)     !
  !                                                                    !
  !   Compressed Row Format:                                           !
  !   - In CRS, a matrix is represented by three arrays: 'val' for     !
  !     non-zero values, 'col' for column indices, and 'row' for       !
  !     row pointers. This format is advantageous for matrix-vector    !
  !     multiplication and other matrix operations.                    !
  !                                                                    !
  !   Example:                                                         !
  !   - Illustrates a 4x4 matrix and its CRS representation.           !
  !     The 'dia' array stores diagonal positions, and 'pos' holds     !
  !     positions in the matrix. 'glo' represents global cell          !
  !     numbering, differing from T-Flows' internal numbering.         !
  !                                                                    !
  !       c   c  .    c                                                !
  !       o   o  .    o                                                !
  !       l   l       l                                                !
  !                                                                    !
  !       1   2       n                                                !
  !                                                                    !
  !     [ 10   0   4   5 ]  --> row 1                                  !
  !     [  2  12  -1   0 ]  --> rows store discretized control volumes !
  !     [  0   1  99   7 ]  ...                                        !
  !     [ -3  11   0  53 ]  --> row n                                  !
  !                                                                    !
  !     Compressed row storage of the above matrix reads:              !
  !                                                                    !
  !     A % val = [  10   4   5   2  12  -1   1  99   7  -3  11  53 ]  !
  !     A % col = [   1   3   4   1   2   3   2   3   4   1   2   4 ]  !
  !     A % row = [   1   4   7  10 ]                                  !
  !                                                                    !
  !     A % dia = [   1   5   9  12 ]                                  !
  !--------------------------------------------------------------------!
  type Sparse_Type
  !> A type which encapsulates matrix in compressed row storage format.

    type(Grid_Type), pointer :: pnt_grid  !! pointer to grid

    integer              :: nonzeros  !! number of nonzero entries
    real,    allocatable :: val(:)    !! value
    real,    allocatable :: fc (:)    !! bare matrix entry for face
    integer, allocatable :: col(:)    !! column positions
    integer, allocatable :: row(:)    !! beginning of each row
    integer, allocatable :: dia(:)    !! diagonal positions
    integer, allocatable :: pos(:,:)  !! face-based position of the matrix
    real,    allocatable :: v_m(:)    !! cell volume over momentum diagonal
    real,    allocatable :: d_inv(:)  !! inverse dia. for preconditioner

    contains
      procedure :: Create_Sparse
      procedure :: Create_Sparse_From_Sparse

  end type

  contains
#   include "Sparse_Mod/Create_Sparse.f90"
#   include "Sparse_Mod/Create_Sparse_From_Sparse.f90"

  end module
