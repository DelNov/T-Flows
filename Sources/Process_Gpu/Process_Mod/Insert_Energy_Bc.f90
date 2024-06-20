!==============================================================================!
  subroutine Insert_Energy_Bc(Process, Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:), cond(:)
  real                      :: a12
  integer                   :: reg, s, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Energy_Bc')

  ! Take some aliases
  b    => Flow % Nat % b
  fc   => Flow % Nat % C % fc
  cond => Flow % conductivity

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$acc parallel loop independent  &
      !$acc present(grid_faces_c, grid_region_f_face, grid_region_l_face,  &
      !$acc         flow_t_n)
      do s = Faces_In_Region_Gpu(reg)
        c1 = grid_faces_c(1,s)
        c2 = grid_faces_c(2,s)
        a12 = cond(c1) * fc(s)
        b(c1) = b(c1) + a12 * flow_t_n(c2)
      end do
      !$acc end parallel

    end if
  end do

  call Profiler % Stop('Insert_Energy_Bc')

  end subroutine
