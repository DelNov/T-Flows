!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(flow)
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
  call Field_Mod_Grad_Variable(flow, flow % u)
  call Field_Mod_Grad_Variable(flow, flow % v)
  call Field_Mod_Grad_Variable(flow, flow % w)

  do c = -grid % n_bnd_cells, grid % n_cells
    flow % shear(c) = u % x(c)**2                     &
                    + v % y(c)**2                     &
                    + w % z(c)**2                     &
                    + 0.5 * (v % z(c) + w % y(c))**2  &
                    + 0.5 * (u % z(c) + w % x(c))**2  &
                    + 0.5 * (v % x(c) + u % y(c))**2

    flow % vort(c) = - (   0.5 * (v % z(c) - w % y(c))**2  &
                         + 0.5 * (u % z(c) - w % x(c))**2  &
                         + 0.5 * (v % x(c) - u % y(c))**2 )

    flow % shear(c) = sqrt(2.0 * flow % shear(c))
    flow % vort(c)  = sqrt(2.0 * abs(flow % vort(c)))
  end do

  end subroutine
