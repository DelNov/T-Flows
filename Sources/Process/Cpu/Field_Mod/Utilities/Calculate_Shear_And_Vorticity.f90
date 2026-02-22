!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(Flow)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the shear stress.                                !
!------------------------------------------------------------------------------!
!   Shear of the velocity vield is a tensor define as:                         !
!                                                                              !
!   s_ij  = 1/2 ( dui/dxj + duj/dxi )                                          !
!                                                                              !
!   Shear's magnitude is computed as:                                          !
!                                                                              !
!   shear = sqrt( 2 * s_ij * s_ij )                                            !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!   Rotation of a velocity field is defined as:                                !
!                                                                              !
!         1         1                                                          !
!   rot = - ∇ x u = - ((dw/dy-dv/dz) i + (du/dz-dw/dx) j + (dv/dx-du/dy) k)    !
!         2         2                                                          !
!                                                                              !
!   Vorticity is twice the rotation vector, hence                              !
!                                                                              !
!                                                                              !
!   vort = ∇ x u = (dw/dy-dv/dz) i + (du/dz-dw/dx) j + (dv/dx-du/dy) k         !
!                                                                              !
!                = v_x i + v_y j + v_z k                                       !
!                                                                              !
!   Vorticity magnitude would be the magnitude of vorticity vector             !
!                                                                              !
!   |vort| = sqrt(v_x^2 + v_y^2 + v_z^2)                                       !
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

  do c = Cells_In_Domain_And_Buffers()
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
