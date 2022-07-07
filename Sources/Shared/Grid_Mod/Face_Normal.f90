!==============================================================================!
  subroutine Face_Normal(Grid, s, nx, ny, nz)
!------------------------------------------------------------------------------!
!   Returns three components of the surface normal vector                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)     :: Grid
  integer, intent(in)  :: s
  real,    intent(out) :: nx
  real,    intent(out) :: ny
  real,    intent(out) :: nz
!==============================================================================!

  nx = Grid % sx(s) / Grid % s(s)
  ny = Grid % sy(s) / Grid % s(s)
  nz = Grid % sz(s) / Grid % s(s)

  end subroutine
