!==============================================================================!
  real function Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump  &
                (mult, phi, s)
!------------------------------------------------------------------------------!
!   Computes cell-face value by gradient interpolation when there is a jump    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type) :: mult
  type(Var_Type)        :: phi
  integer               :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  real                      :: w1, w2, dx1, dy1, dz1, dx2, dy2, dz2
  integer                   :: c1, c2
!==============================================================================!

  flow => mult % pnt_flow
  grid => flow % pnt_grid

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  ! Vector from c1 to face
  dx1 = grid % xf(s) - grid % xc(c1)
  dy1 = grid % yf(s) - grid % yc(c1)
  dz1 = grid % zf(s) - grid % zc(c1)

  ! Vector from c2 to face
  dx2 = grid % xf(s) - (grid % xc(c1) + grid % dx(s))
  dy2 = grid % yf(s) - (grid % yc(c1) + grid % dy(s))
  dz2 = grid % zf(s) - (grid % zc(c1) + grid % dz(s))

  w1 = grid % f(s)
  if(mult % cell_at_elem(c1) .eq. 0 .and.  &
     mult % cell_at_elem(c2) .ne. 0) w1 = 1.0
  if(mult % cell_at_elem(c2) .eq. 0 .and.  &
     mult % cell_at_elem(c1) .ne. 0) w1 = 0.0
  w2 = 1.0 - w1

  Multiphase_Mod_Vof_Interpolate_Var_To_Face_With_Jump                        &
                                    = w1 * (phi % n(c1) + phi % x(c1) * dx1   &
                                                        + phi % y(c1) * dy1   &
                                                        + phi % z(c1) * dz1)  &
                                    + w2 * (phi % n(c2) + phi % x(c2) * dx2   &
                                                        + phi % y(c2) * dy2   &
                                                        + phi % z(c2) * dz2)

  end function

