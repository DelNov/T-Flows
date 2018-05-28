!==============================================================================!
  subroutine Allocate_Memory(grid)
!------------------------------------------------------------------------------!
! Alocates memory for geometrical quantities.                                  !
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Comm_Mod
  use Grid_Mod
  use Work_Mod
  use Matrix_Mod
  use Solvers_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!==============================================================================!

  ! Variables defined in all.h90:
  allocate (fw(grid % n_faces));  fw = 0.0  

  ! Variables defined in Solvers_Mod
  call Matrix_Mod_Allocate(grid, D)

  ! Variables defined in Flow_Mod.h90:
  call Matrix_Mod_Allocate(grid, a)
  allocate (b(grid % n_cells));  b=0.

  ! Working arrays
  call Work_Mod_Allocate_Real_Cells(grid, 30)
  call Work_Mod_Allocate_Real_Faces(grid,  1)
  call Work_Mod_Allocate_Real_Nodes(grid,  1)

  call Work_Mod_Allocate_Integer_Cells(grid, 4)

  ! Variables defined in Comm_Mod.h90:
  call Comm_Mod_Allocate_Memory(grid)

  ! This array should be handled in a more elegant way
  allocate (f_coef(grid % n_faces)); f_coef=0.

  end subroutine
