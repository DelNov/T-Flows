!==============================================================================!
  subroutine Add_Body_Force_Term(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, contiguous, pointer :: ui_o(:)
  real, contiguous, pointer :: b(:), vol(:)
  integer                   :: c, nc
  real                      :: drop
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Body_Force_Term')

  ! Take some aliases
  Grid => Flow % pnt_grid
  b    => Flow % Nat % b
  vol  => Grid % vol
  nc   =  Grid % n_cells

  if(comp .eq. 1) ui_o => Flow % u % o
  if(comp .eq. 2) ui_o => Flow % v % o
  if(comp .eq. 3) ui_o => Flow % w % o
  if(comp .eq. 1) drop =  Flow % bulk % p_drop_x
  if(comp .eq. 2) drop =  Flow % bulk % p_drop_y
  if(comp .eq. 3) drop =  Flow % bulk % p_drop_z

  !$acc parallel loop independent
  do c = 1, nc
    b(c) = b(c) + drop * vol(c)
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Body_Force_Term')

  end subroutine
