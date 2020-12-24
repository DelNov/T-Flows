!==============================================================================!
  subroutine Grid_Mod_Bounding_Box(grid, xmin, ymin, zmin, xmax, ymax, zmax)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: xmin, ymin, zmin, xmax, ymax, zmax
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn
!==============================================================================!

  ! Take alias
  nn = grid % n_nodes

  xmin = minval(grid % xn(1:nn))
  ymin = minval(grid % yn(1:nn))
  zmin = minval(grid % zn(1:nn))

  xmax = maxval(grid % xn(1:nn))
  ymax = maxval(grid % yn(1:nn))
  zmax = maxval(grid % zn(1:nn))

  call Comm_Mod_Global_Min_Real(xmin)
  call Comm_Mod_Global_Min_Real(ymin)
  call Comm_Mod_Global_Min_Real(zmin)

  call Comm_Mod_Global_Max_Real(xmax)
  call Comm_Mod_Global_Max_Real(ymax)
  call Comm_Mod_Global_Max_Real(zmax)

  end subroutine
