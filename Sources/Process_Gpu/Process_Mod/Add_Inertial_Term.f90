!==============================================================================!
  subroutine Add_Inertial_Term(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: ui_o(:), b(:)
  integer                   :: c
  real                      :: dens, dt
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  b    => Flow % Nat % b
  dens =  Flow % density
  dt   =  Flow % dt

  if(comp .eq. 1) ui_o => Flow % u % o
  if(comp .eq. 2) ui_o => Flow % v % o
  if(comp .eq. 3) ui_o => Flow % w % o

  !$acc parallel loop independent
  do c = 1, grid_n_cells
    b(c) = b(c) + dens * ui_o(c) * grid_vol(c) / dt
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
