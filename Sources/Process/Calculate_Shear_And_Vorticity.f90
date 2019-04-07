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
  use Grad_Mod
  use Field_Mod,  only: Field_Type
  use Grid_Mod,   only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  ! Velocity gradients
  call Grad_Mod_Variable(flow % u, .true.)
  call Grad_Mod_Variable(flow % v, .true.)
  call Grad_Mod_Variable(flow % w, .true.)

  shear(:) =  u % x(:)**2                     &
            + v % y(:)**2                     &
            + w % z(:)**2                     &
            + 0.5 * (v % z(:) + w % y(:))**2  &
            + 0.5 * (u % z(:) + w % x(:))**2  &
            + 0.5 * (v % x(:) + u % y(:))**2

  vort(:) = - (   0.5*(v % z(:) - w % y(:))**2  &
                + 0.5*(u % z(:) - w % x(:))**2  &
                + 0.5*(v % x(:) - u % y(:))**2 )

  shear(:) = sqrt(2.0 * shear(:))
  vort(:)  = sqrt(2.0 * abs(vort(:)))

  end subroutine
