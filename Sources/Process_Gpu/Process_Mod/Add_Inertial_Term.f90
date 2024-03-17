!==============================================================================!
  subroutine Add_Inertial_Term(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: ui_o(:), b(:), dens(:)
  integer                   :: c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  b    => Flow % Nat % b
  dens => Flow % density

  if(comp .eq. 1) ui_o => Flow % u % o
  if(comp .eq. 2) ui_o => Flow % v % o
  if(comp .eq. 3) ui_o => Flow % w % o

  !$acc parallel loop independent
  do c = 1, grid_n_cells - grid_n_buff_cells
    b(c) = b(c) + dens(c) * ui_o(c) * grid_vol(c) / Flow % dt
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
