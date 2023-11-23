#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Matrix_Mod
!------------------------------------------------------------------------------!
!>  The Matrix_Mod module is designed to handle the storage and manipulation of
!>  sparse matrices within T-Flows. It stores matrices in the compressed row
!>  storage (CSR) format, a common format for storing sparse matrices.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !   Matrix type                                                      !
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
  type Matrix_Type
  !> A type which encapsulates matrix in compressed row storage format.

    type(Grid_Type), pointer :: pnt_grid  !! pointer to grid

    integer              :: nonzeros  !! number of nonzero entries
    real,    allocatable :: val(:)    !! value
    real,    allocatable :: fc (:)    !! bare matrix entry for face
    real,    allocatable :: sav(:)    !! saved momentum diagonal value
    integer, allocatable :: col(:)    !! beginning of each row
    integer, allocatable :: row(:)    !! column positions
    integer, allocatable :: dia(:)    !! diagonal positions
    integer, allocatable :: pos(:,:)  !! position in the matrix
    integer, allocatable :: glo(:)    !! global cell number for PETSc
                                      !! which is different from T-Flows
                                      !! and stars from zero
    contains
      procedure :: Create_Matrix

  end type

  contains

#   include "Matrix_Mod/Create_Matrix.f90"

  end module
