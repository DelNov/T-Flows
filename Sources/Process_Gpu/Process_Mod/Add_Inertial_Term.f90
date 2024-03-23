!==============================================================================!
  subroutine Add_Inertial_Term(Process, Flow, Grid, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: ui_o(:), b(:), dens(:)
  integer                   :: c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  b    => Flow % Nat % b
  dens => Flow % density

  if(comp .eq. 1) ui_o => Flow % u % o
  if(comp .eq. 2) ui_o => Flow % v % o
  if(comp .eq. 3) ui_o => Flow % w % o

  !$acc parallel loop independent
  do c = Cells_In_Domain()
    b(c) = b(c) + dens(c) * ui_o(c) * Grid % vol(c) / Flow % dt
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
