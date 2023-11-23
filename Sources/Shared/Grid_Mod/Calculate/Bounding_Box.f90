!==============================================================================!
  subroutine Bounding_Box(Grid, xmin, ymin, zmin, xmax, ymax, zmax)
!------------------------------------------------------------------------------!
!>  Returns the extens of a grid's bounding box.  It is called only from
!>  Print_Grid_Statistics which, in turn, is called only from Convert.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid                                !! parent grid class
  real             :: xmin, ymin, zmin, xmax, ymax, zmax  !! coordinate extent
!-----------------------------------[Locals]-----------------------------------!
  integer :: nn
!==============================================================================!

  ! Take alias
  nn = Grid % n_nodes

  xmin = minval(Grid % xn(1:nn))
  ymin = minval(Grid % yn(1:nn))
  zmin = minval(Grid % zn(1:nn))

  xmax = maxval(Grid % xn(1:nn))
  ymax = maxval(Grid % yn(1:nn))
  zmax = maxval(Grid % zn(1:nn))

  call Global % Min_Real(xmin)
  call Global % Min_Real(ymin)
  call Global % Min_Real(zmin)

  call Global % Max_Real(xmax)
  call Global % Max_Real(ymax)
  call Global % Max_Real(zmax)

  end subroutine
