!==============================================================================!
  subroutine Calculate_Vorticity(flow)
!------------------------------------------------------------------------------!
!  Computes the magnitude of the vorticity                                     !
!------------------------------------------------------------------------------!
!  vort = sqrt( 2 * Sij * Sij )                                                !
!  Sij = 1/2 ( dUi/dXj - dUj/dXi )                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
  use Comm_Mod
  use Turb_Mod
  use Grid_Mod
  use Grad_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! Velocity gradients
  call Grad_Mod_Variable(flow % u)
  call Grad_Mod_Variable(flow % v)
  call Grad_Mod_Variable(flow % w)

  do c = 1, grid % n_cells
    flow % vort(c) = 2.0 * (0.5 * (w % y(c) - v % z(c)))**2  &
                   + 2.0 * (0.5 * (w % x(c) - u % z(c)))**2  &
                   + 2.0 * (0.5 * (v % x(c) - u % y(c)))**2

    flow % vort(c) = sqrt(abs(2.0 * flow % vort(c)))
  end do

  end subroutine
