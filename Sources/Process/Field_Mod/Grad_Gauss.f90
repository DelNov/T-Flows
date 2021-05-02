!==============================================================================!
  subroutine Field_Mod_Grad_Gauss(flow, phi_f, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradients in cells from values in face with Gauss' theorem.     !
!   It heavily relies on the accuracy of values in faces, which are not        !
!   calculated here, and the primary use of this function is to be used as     !
!   embedded in an iterative algorithm which also updates face values.         !
!   See also it's sister function Interpolate_To_Faces from this module, and   !
!   its parent function Grad_Gauss_Variable, also from this module.            !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phi_f(1:flow % pnt_grid % n_faces)
  real                     :: phi_x( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
  real                     :: phi_y( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
  real                     :: phi_z( -flow % pnt_grid % n_bnd_cells:  &
                                      flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2, c
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  !-----------------------------------------------!
  !   Update gradients from the values at faces   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
  end do

  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then
      phi_x(c1) = phi_x(c1) - phi_f(s) * grid % sx(s)
      phi_y(c1) = phi_y(c1) - phi_f(s) * grid % sy(s)
      phi_z(c1) = phi_z(c1) - phi_f(s) * grid % sz(s)

      phi_x(c2) = phi_x(c2) + phi_f(s) * grid % sx(s)
      phi_y(c2) = phi_y(c2) + phi_f(s) * grid % sy(s)
      phi_z(c2) = phi_z(c2) + phi_f(s) * grid % sz(s)
    else
      phi_x(c1) = phi_x(c1) - phi_f(s) * grid % sx(s)
      phi_y(c1) = phi_y(c1) - phi_f(s) * grid % sy(s)
      phi_z(c1) = phi_z(c1) - phi_f(s) * grid % sz(s)
    end if

  end do

  do c = 1, grid % n_cells
    phi_x(c) = -phi_x(c) / grid % vol(c)
    phi_y(c) = -phi_y(c) / grid % vol(c)
    phi_z(c) = -phi_z(c) / grid % vol(c)
  end do

  end subroutine
