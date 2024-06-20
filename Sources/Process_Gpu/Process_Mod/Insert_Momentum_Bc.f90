!==============================================================================!
  subroutine Insert_Momentum_Bc(Process, Flow, Grid, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type)          :: Grid
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:), ui_n(:), visc(:)
  real                      :: m12
  integer                   :: reg, s, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Momentum_Bc')

  ! Take some aliases
  b    => Flow % Nat % b
  fc   => Flow % Nat % C % fc
  visc => Flow % viscosity

  if(comp .eq. 1) ui_n => flow_u_n
  if(comp .eq. 2) ui_n => flow_v_n
  if(comp .eq. 3) ui_n => flow_w_n

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$acc parallel loop independent                        &
      !$acc present(grid_region_f_face, grid_region_l_face,  &
      !$acc         grid_faces_c,                            &
      !$acc         visc, fc, b, ui_n)
      do s = Faces_In_Region_Gpu(reg)  ! all present
        c1 = grid_faces_c(1,s)
        c2 = grid_faces_c(2,s)
        m12 = visc(c1) * fc(s)
        b(c1) = b(c1) + m12 * ui_n(c2)
      end do
      !$acc end parallel

    end if
  end do

  call Profiler % Stop('Insert_Momentum_Bc')

  end subroutine
