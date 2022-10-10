!==============================================================================!
  subroutine Calculate_Vorticity(Flow)
!------------------------------------------------------------------------------!
!  Computes the magnitude of the vorticity                                     !
!------------------------------------------------------------------------------!
!  vort = sqrt( 2 * Sij * Sij )                                                !
!  Sij = 1/2 ( dUi/dXj - dUj/dXi )                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Velocity gradients
  call Flow % Grad_Variable(Flow % u)
  call Flow % Grad_Variable(Flow % v)
  call Flow % Grad_Variable(Flow % w)

  do c = 1, Grid % n_cells
    Flow % vort(c) = 2.0 * (0.5 * (w % y(c) - v % z(c)))**2  &
                   + 2.0 * (0.5 * (w % x(c) - u % z(c)))**2  &
                   + 2.0 * (0.5 * (v % x(c) - u % y(c)))**2

    Flow % vort(c) = sqrt(abs(2.0 * Flow % vort(c)))
  end do

  end subroutine
