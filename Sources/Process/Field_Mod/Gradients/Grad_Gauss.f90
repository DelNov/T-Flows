!==============================================================================!
  subroutine Grad_Gauss(Flow, phi_f, phi_x, phi_y, phi_z)
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
  class(Field_Type), target :: Flow
  real                      :: phi_f(1:Flow % pnt_grid % n_faces)
  real                      :: phi_x( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
  real                      :: phi_y( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
  real                      :: phi_z( -Flow % pnt_grid % n_bnd_cells:  &
                                       Flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2, c, reg
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  !-----------------------------------------------!
  !   Update gradients from the values at faces   !
  !-----------------------------------------------!
  do c = Cells_In_Domain()
    phi_x(c) = 0.0
    phi_y(c) = 0.0
    phi_z(c) = 0.0
  end do

  ! Faces in boundary region first
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)

      phi_x(c1) = phi_x(c1) - phi_f(s) * Grid % sx(s)
      phi_y(c1) = phi_y(c1) - phi_f(s) * Grid % sy(s)
      phi_z(c1) = phi_z(c1) - phi_f(s) * Grid % sz(s)
    end do
  end do

  ! Faces inside the domain
  do s = Faces_In_Domain()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_x(c1) = phi_x(c1) - phi_f(s) * Grid % sx(s)
    phi_y(c1) = phi_y(c1) - phi_f(s) * Grid % sy(s)
    phi_z(c1) = phi_z(c1) - phi_f(s) * Grid % sz(s)

    phi_x(c2) = phi_x(c2) + phi_f(s) * Grid % sx(s)
    phi_y(c2) = phi_y(c2) + phi_f(s) * Grid % sy(s)
    phi_z(c2) = phi_z(c2) + phi_f(s) * Grid % sz(s)
  end do

  do c = Cells_In_Domain()
    phi_x(c) = -phi_x(c) / Grid % vol(c)
    phi_y(c) = -phi_y(c) / Grid % vol(c)
    phi_z(c) = -phi_z(c) / Grid % vol(c)
  end do

  end subroutine
