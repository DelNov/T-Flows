!==============================================================================!
  real function Interpolate_Var_To_Face_With_Jump(Vof, phi, s)
!------------------------------------------------------------------------------!
!   Computes cell-face value by gradient interpolation when there is a jump    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type) :: Vof
  type(Var_Type)        :: phi
  integer               :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  real                      :: w1, w2, dx1, dy1, dz1, dx2, dy2, dz2
  integer                   :: c1, c2
!==============================================================================!

  Flow => Vof % pnt_flow
  Grid => Flow % pnt_grid

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  ! Vector from c1 to face
  dx1 = Grid % xf(s) - Grid % xc(c1)
  dy1 = Grid % yf(s) - Grid % yc(c1)
  dz1 = Grid % zf(s) - Grid % zc(c1)

  ! Vector from c2 to face
  dx2 = Grid % xf(s) - (Grid % xc(c1) + Grid % dx(s))
  dy2 = Grid % yf(s) - (Grid % yc(c1) + Grid % dy(s))
  dz2 = Grid % zf(s) - (Grid % zc(c1) + Grid % dz(s))

  w1 = Grid % f(s)
  if(Vof % cell_at_elem(c1) .eq. 0 .and.  &
     Vof % cell_at_elem(c2) .ne. 0) w1 = 1.0
  if(Vof % cell_at_elem(c2) .eq. 0 .and.  &
     Vof % cell_at_elem(c1) .ne. 0) w1 = 0.0
  w2 = 1.0 - w1

  Interpolate_Var_To_Face_With_Jump                        &
                 = w1 * (phi % n(c1) + phi % x(c1) * dx1   &
                                     + phi % y(c1) * dy1   &
                                     + phi % z(c1) * dz1)  &
                 + w2 * (phi % n(c2) + phi % x(c2) * dx2   &
                                     + phi % y(c2) * dy2   &
                                     + phi % z(c2) * dz2)

  end function

