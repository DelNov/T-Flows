!==============================================================================!
  subroutine Face_Normal(Grid, s, nx, ny, nz)
!------------------------------------------------------------------------------!
!>  Returns three components of the surface normal vector.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)     :: Grid  !! grid under consideration
  integer, intent(in)  :: s     !! face number
  real,    intent(out) :: nx    !! x component of the surface normal
  real,    intent(out) :: ny    !! y component of the surface normal
  real,    intent(out) :: nz    !! z compomnent of the surface normal
!==============================================================================!

  nx = Grid % sx(s) / Grid % s(s)
  ny = Grid % sy(s) / Grid % s(s)
  nz = Grid % sz(s) / Grid % s(s)

  end subroutine
