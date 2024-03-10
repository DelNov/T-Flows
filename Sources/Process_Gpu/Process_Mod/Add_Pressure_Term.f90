!==============================================================================!
  subroutine Add_Pressure_Term(Proc, Flow, comp)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Dimension of the system under consideration                                !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!                                                                              !
!   Pressure gradient alone:                                                   !
!     p % x        [kg/(m^2 s^2)]                                              !
!                                                                              !
!   Pressure gradient times volume:                                            !
!     p % x * vol  [kg/(m^2 s^2) * m^3 = kg m / s^2 = N]                       !
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, contiguous, pointer :: b(:), p_i(:), vol(:)
  real                      :: p_d_i
  integer                   :: c, nc
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Pressure_Term')

  ! Take some aliases
  Grid => Flow % pnt_grid
  b    => Flow % Nat % b
  vol  => Grid % vol
  nc   =  Grid % n_cells

  ! Still on aliases
  if(comp .eq. 1) p_i   => Flow % p % x
  if(comp .eq. 2) p_i   => Flow % p % y
  if(comp .eq. 3) p_i   => Flow % p % z
  if(comp .eq. 1) p_d_i =  Flow % p_drop_x
  if(comp .eq. 2) p_d_i =  Flow % p_drop_y
  if(comp .eq. 3) p_d_i =  Flow % p_drop_z

  !$acc parallel loop independent
  do c = 1, nc
    b(c) = b(c) + (p_d_i - p_i(c)) * vol(c)
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Pressure_Term')

  end subroutine
