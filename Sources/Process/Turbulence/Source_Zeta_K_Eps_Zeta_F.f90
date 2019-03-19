!==============================================================================!
  subroutine Source_Zeta_K_Eps_Zeta_F(flow, sol, n_step)
!------------------------------------------------------------------------------!
!   Calculate source terms in equation for zeta.                               !
!   Term which is negative is put on left hand side in diagonal of             !
!   martix of coefficient                                                      !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Field_Mod
  use Rans_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  integer                   :: n_step
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c
!==============================================================================!
!   In transport equation for zeta two source terms exist which have form:     !
!                                                                              !
!    /                                                                         !
!    |                     eps                                                 !
!    | ( f22 * kin * dV - ---- * dV )                                          !
!    |                    zeta                                                 !
!    /                                                                         !
!                                                                              !
!   First term can appear as positive and as negative as well so depend of     !
!   sign of term , it is placed on left or right hand side.  Second, negative  !
!   source term is added to main diagonal left hand side coefficient matrix    !
!   in order to increase stability of solver                                   !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  a    => sol % a
  b    => sol % b % val

  ! Positive source term 
  ! The first option in treating the source is making computation very 
  ! sensitive to initial condition while the second one can lead to 
  ! instabilities for some cases such as flow around cylinder. That is why we
  ! choose this particular way to the add source term.
  do c = 1, grid % n_cells
    if(n_step > 500) then
      b(c) = b(c) + f22 % n(c) * grid % vol(c) * density
    else
      b(c) = b(c) + max(0.0, f22 % n(c)*grid % vol(c)) * density
      a % val(a % dia(c)) = a % val(a % dia(c))                  &
                          + max(0.0, -f22 % n(c) * grid % vol(c) &
                          / (zeta % n(c) + TINY)) * density
    end if
    a % val(a % dia(c)) =  a % val(a % dia(c))      &
                        + grid % vol(c) * p_kin(c)  &
                        / (kin % n(c) + TINY) 
  end do

  end subroutine
