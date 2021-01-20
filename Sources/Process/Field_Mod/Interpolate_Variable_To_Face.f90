!==============================================================================!
  real function Field_Mod_Interpolate_Variable_To_Face(flow, phi, s)
!------------------------------------------------------------------------------!
!   Computes cell-face value by gradient interpolation.                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: phi
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

  Field_Mod_Interpolate_Variable_To_Face =                 &
      phi % n(c1) * w1 + phi % n(c2) * w2                  &
   + (phi % x(c1) * w1 + phi % x(c2) * w2) * grid % rx(s)  &
   + (phi % y(c1) * w1 + phi % y(c2) * w2) * grid % ry(s)  &
   + (phi % z(c1) * w1 + phi % z(c2) * w2) * grid % rz(s)

  end function

