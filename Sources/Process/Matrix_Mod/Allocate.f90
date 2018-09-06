!==============================================================================!
  subroutine Matrix_Mod_Allocate(grid, matrix)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)   :: grid
  type(Matrix_Type) :: matrix
!==============================================================================!

  ! Allocate memory for matrix
  allocate (matrix % row(  grid % n_cells + 1));                matrix % row=0
  allocate (matrix % dia(  grid % n_cells));                    matrix % dia=0
  allocate (matrix % sav( -grid % n_bnd_cells:grid % n_cells)); matrix % sav=0.
  allocate (matrix % pos(2,grid % n_faces));                    matrix % pos=0

  end subroutine
