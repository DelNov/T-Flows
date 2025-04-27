!==============================================================================!
  subroutine Create_Sparse_Val(Aval, Acon)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Val_Type) :: Aval  !! parent class
  type (Sparse_Con_Type) :: Acon  !! connectivity matrix
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Fetch the alias
  Grid => Acon % pnt_grid

  ! Store pointer to the grid
  Aval % pnt_grid => Grid

  O_Print '(a)', ' # Creating a sparse value matrix'

  allocate(Aval % val  (Acon % nonzeros));  Aval % val   = 0
  allocate(Aval % d_inv(Grid % n_cells));   Aval % d_inv = 0
  
  end subroutine
