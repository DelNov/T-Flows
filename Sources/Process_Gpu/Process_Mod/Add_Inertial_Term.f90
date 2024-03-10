!==============================================================================!
  subroutine Add_Inertial_Term(Proc, Flow, comp)
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
  real                      :: dens, dt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  Grid => Flow % pnt_grid
  b    => Flow % Nat % b
  vol  => Grid % vol
  nc   =  Grid % n_cells
  dens =  Flow % density
  dt   =  Flow % dt

  if(comp .eq. 1) ui_o => Flow % u % o
  if(comp .eq. 2) ui_o => Flow % v % o
  if(comp .eq. 3) ui_o => Flow % w % o

  !$acc parallel loop independent
  do c = 1, nc
    b(c) = b(c) + dens * ui_o(c) * vol(c) / dt
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
