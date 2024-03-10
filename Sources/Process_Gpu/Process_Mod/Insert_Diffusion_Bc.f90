!==============================================================================!
  subroutine Insert_Diffusion_Bc(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  real,            pointer :: b(:), ui(:)
  real                     :: visc, m12
  integer                  :: reg, s, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Insert_Diffusion_Bc')

  ! Take some aliases
  Grid => Flow % pnt_grid
  b    => Flow % Nat % b
  visc =  Flow % viscosity
  if(comp .eq. 1) ui => Flow % u % n
  if(comp .eq. 2) ui => Flow % v % n
  if(comp .eq. 3) ui => Flow % w % n

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL) then
      !$acc parallel loop independent
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        m12 = visc * Grid % s(s) / Grid % d(s)

        b(c1) = b(c1) + m12 * ui(c2)
      end do
      !$acc end parallel
    end if
  end do

  call Profiler % Stop('Insert_Diffusion_Bc')

  end subroutine
