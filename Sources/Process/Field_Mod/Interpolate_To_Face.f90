!==============================================================================!
  real function Field_Mod_Interpolate_To_Face(flow, phi, phi_x, phi_y, phi_z, s)
!------------------------------------------------------------------------------!
!   Computes cell-face value by gradient interpolation.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  real             :: phi  ( -flow % pnt_grid % n_bnd_cells  &
                             :flow % pnt_grid % n_cells)
  real             :: phi_x( -flow % pnt_grid % n_bnd_cells  &
                             :flow % pnt_grid % n_cells)
  real             :: phi_y( -flow % pnt_grid % n_bnd_cells  &
                             :flow % pnt_grid % n_cells)
  real             :: phi_z( -flow % pnt_grid % n_bnd_cells  &
                             :flow % pnt_grid % n_cells)
  integer          :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  real                     :: w1, w2
  integer                  :: c1, c2
!==============================================================================!

  grid => flow % pnt_grid

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  w1 = grid % f(s)
  w2 = 1.0 - w1

  Field_Mod_Interpolate_To_Face =                      &
      phi(c1)   * w1 + phi(c2)   * w2                  &
   + (phi_x(c1) * w1 + phi_x(c2) * w2) * grid % rx(s)  &
   + (phi_y(c1) * w1 + phi_y(c2) * w2) * grid % ry(s)  &
   + (phi_z(c1) * w1 + phi_z(c2) * w2) * grid % rz(s)

  end function

