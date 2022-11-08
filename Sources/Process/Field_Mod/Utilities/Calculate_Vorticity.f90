!==============================================================================!
  subroutine Calculate_Vorticity(Flow)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the vorticity                                    !
!------------------------------------------------------------------------------!
!   Rotation of a velocity field is defined as:                                !
!                                                                              !
!         1         1     dw   dv         du   dw         dv   du              !
!   rot = - ∇ x u = - ( ( -- - -- ) i + ( -- - -- ) j + ( -- - -- ) k )        !
!         2         2     dy   dz         dz   dx         dx   dz              !
!                                                                              !
!   Vorticity is twice the rotation vector, hence                              !
!                                                                              !
!                    dw   dv         du   dw         dv   du                   !
!   vort = ∇ x u = ( -- - -- ) i + ( -- - -- ) j + ( -- - -- ) k )             !
!                    dy   dz         dz   dx         dx   dz                   !
!                                                                              !
!                = v_x i + v_y j + v_z k                                       !
!                                                                              !
!   Vorticity magnitude would be the magnitude of vorticity vector             !
!                                                                              !
!   |vort| = sqrt(v_x^2 + v_y^2 + v_z^2)                                       !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
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
