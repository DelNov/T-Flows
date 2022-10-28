!==============================================================================!
  subroutine Allocate_Work(Work, Grid, n_r_cell, n_r_face, n_r_node,  &
                                       n_i_cell, n_i_face, n_i_node)
!------------------------------------------------------------------------------!
!   Alocates memory for working arrays                                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work
  type(Grid_Type)  :: Grid(:)
  integer          :: n_r_cell  ! number of real cell arrays
  integer          :: n_r_face  ! number of real face arrays
  integer          :: n_r_node  ! number of real node arrays
  integer          :: n_i_cell  ! number of integer cell arrays
  integer          :: n_i_face  ! number of integer face arrays
  integer          :: n_i_node  ! number of integer node arrays
!==============================================================================!

  Work % req_r_cell = n_r_cell
  Work % req_r_face = n_r_face
  Work % req_r_node = n_r_node

  Work % req_i_cell = n_i_cell
  Work % req_i_face = n_i_face
  Work % req_i_node = n_i_node

  call Work % Allocate_Real_Cell(Grid, n_r_cell)
  call Work % Allocate_Real_Face(Grid, n_r_face)
  call Work % Allocate_Real_Node(Grid, n_r_node)

  call Work % Allocate_Int_Cell(Grid, n_i_cell)
  call Work % Allocate_Int_Face(Grid, n_i_face)
  call Work % Allocate_Int_Node(Grid, n_i_node)

  end subroutine
