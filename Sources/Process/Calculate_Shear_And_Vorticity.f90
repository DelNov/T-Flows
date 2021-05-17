!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(Flow)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the shear stress.                                !
!------------------------------------------------------------------------------!
!   shear = sqrt( 2 * Sij * Sij )                                              !
!   Sij = 1/2 ( dUi/dXj + dUj/dXi )                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Turb_Mod
  use Field_Mod,  only: Field_Type
  use Grid_Mod,   only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  ! Velocity gradients
  call Flow % Grad_Variable(Flow % u)
  call Flow % Grad_Variable(Flow % v)
  call Flow % Grad_Variable(Flow % w)

  do c = -grid % n_bnd_cells, grid % n_cells
    Flow % shear(c) = u % x(c)**2                     &
                    + v % y(c)**2                     &
                    + w % z(c)**2                     &
                    + 0.5 * (v % z(c) + w % y(c))**2  &
                    + 0.5 * (u % z(c) + w % x(c))**2  &
                    + 0.5 * (v % x(c) + u % y(c))**2

    Flow % vort(c) = - (   0.5 * (v % z(c) - w % y(c))**2  &
                         + 0.5 * (u % z(c) - w % x(c))**2  &
                         + 0.5 * (v % x(c) - u % y(c))**2 )

    Flow % shear(c) = sqrt(2.0 * Flow % shear(c))
    Flow % vort(c)  = sqrt(2.0 * abs(Flow % vort(c)))
  end do

  end subroutine
