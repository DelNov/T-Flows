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
  real,            pointer :: b(:), ui_n(:), grid_s(:), grid_d(:)
  integer,         pointer :: grid_faces_c(:,:)
  integer,         pointer :: grid_reg_f_face(:), grid_reg_l_face(:)
  real                     :: visc, m12
  integer                  :: reg, s, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Insert_Diffusion_Bc')

  ! Take some aliases
  Grid             => Flow % pnt_grid
  b                => Flow % Nat % b
  grid_s           => Grid % s
  grid_d           => Grid % d
  grid_reg_f_face  => Grid % region % f_face
  grid_reg_l_face  => Grid % region % l_face
  grid_faces_c     => Grid % faces_c
  visc             =  Flow % viscosity

  if(comp .eq. 1) ui_n => Flow % u % n
  if(comp .eq. 2) ui_n => Flow % v % n
  if(comp .eq. 3) ui_n => Flow % w % n

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!
  !$acc kernels
  b(:) = 0.0
  !$acc end kernels

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL) then

      !$acc parallel loop independent
      do s = grid_reg_f_face(reg), grid_reg_l_face(reg)
        c1 = grid_faces_c(1,s)
        c2 = grid_faces_c(2,s)
        m12 = visc * grid_s(s) / grid_d(s)
        b(c1) = b(c1) + m12 * ui_n(c2)
      end do
      !$acc end parallel

    end if
  end do

  call Profiler % Stop('Insert_Diffusion_Bc')

  end subroutine
